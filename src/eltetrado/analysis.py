import copy
import itertools
import logging
import math
import os
import subprocess
import sys
import tempfile
from collections import Counter, defaultdict
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import IO, TYPE_CHECKING

import numpy
import numpy.typing
from rnapolis.common import (
    BaseInteractions,
    BpSeq,
    Entry,
    GlycosidicBond,
    LeontisWesthof,
    Residue,
)
from rnapolis.tertiary import (
    Atom,
    BasePair3D,
    Mapping2D3D,
    Residue3D,
    Structure3D,
    torsion_angle,
)

from eltetrado.model import (
    ONZ,
    ONZM,
    Direction,
    GbaQuadruplexClassification,
    GbaTetradClassification,
    HelixHandedness,
    Ion,
    LoopClassification,
    LoopType,
    SnapbackOrientation,
    StrandPolarity,
    SugarPucker,
    TetradPolarity,
)

if TYPE_CHECKING:
    from rnapolis.common import MissingResidue

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))
logger = logging.getLogger(__name__)

# Limits of the spurious quadruplex filter.
#
# Tetrad discovery is permissive on purpose: a tetrad only needs a cycle of
# Leontis/Westhof base-pair edges (see ``Tetrad.is_valid``), and two tetrads
# become a stack when as little as a single chain-adjacent or stacked
# nucleotide supports the pairing. In large structured RNAs (e.g. ribosomal
# RNA) this occasionally produces pseudo-quadruplexes: cycles of mixed bases
# with unrealistic stacking geometry, joined by chain segments hundreds of
# nucleotides long. The limits below reject such structures; they were tuned
# on the ElTetrado/onquadro corpus (DSSR-annotated), where real quadruplexes
# have bulges of up to 15 nucleotides (7cls) and pairs with rise between 2.69
# and 10.4 Å (7csk, 8hb3), so keep the limits generous. Twist is deliberately
# NOT limited: real two-block structures (7d5d, 7d5e, 7d5f, 6qjo) contain
# pairs with twist up to 141°, and across antiparallel/snapback steps (2la5,
# 5dea) it reaches 170°, so twist does not separate real stacks from pseudo
# ones. Handedness is not a discriminator either — real left-handed
# quadruplexes exist both as two-block and plain structures.
MIN_QUADRUPLEX_TETRADS = 2
MAX_QUADRUPLEX_BULGES = 25
MIN_QUADRUPLEX_RISE = 2.0
MAX_QUADRUPLEX_RISE = 12.0


def get_plane_vectors(coords: numpy.ndarray) -> tuple[numpy.ndarray, numpy.ndarray]:
    """
    Helper function to get centroid and normal for a single set of coordinates.

    Parameters:
    coords (numpy.ndarray): An N x 3 array of XYZ coordinates.

    Returns:
    tuple: (centroid, normal_vector) where both are numpy arrays
    """
    centroid = numpy.mean(coords, axis=0)
    centered = coords - centroid
    _u, _s, vt = numpy.linalg.svd(centered)
    # The normal is the last row of Vt (smallest variance)
    normal = vt[2, :]
    return centroid, normal


def calculate_planarity_deviation(atom_coords: numpy.ndarray) -> dict[str, float]:
    """
    Calculates the RMSD of atom coordinates from their best-fit plane.

    Parameters:
    atom_coords (numpy.ndarray): An N x 3 array of XYZ coordinates.

    Returns:
    dict: Dictionary with keys 'rmsd', 'max', 'avg' containing:
        - 'rmsd': Root Mean Square Deviation from the plane
        - 'max': Maximum absolute deviation from the plane
        - 'avg': Average absolute deviation from the plane
    """
    centroid, normal_vector = get_plane_vectors(atom_coords)

    # Center the coordinates
    centered_coords = atom_coords - centroid

    # Calculate Distances
    # Project the centered coordinates onto the normal vector.
    # The dot product gives the signed distance from the plane for each atom.
    # d = (P - C) . n
    distances = numpy.dot(centered_coords, normal_vector)

    # Compute RMSD (Root Mean Square Deviation)
    rmsd = numpy.sqrt(numpy.mean(distances**2))
    max_deviation = numpy.max(numpy.abs(distances))
    avg_deviation = numpy.mean(numpy.abs(distances))

    return {"rmsd": rmsd, "max": max_deviation, "avg": avg_deviation}


def collect_nucleobase_atoms(nucleotides: tuple[Residue3D, ...]) -> list[numpy.ndarray]:
    """
    Collect all heavy atoms from nucleobases of given nucleotides.

    Parameters:
    nucleotides (tuple): Tuple of Residue3D objects.

    Returns:
    list: List of coordinate arrays for all nucleobase heavy atoms.
    """
    coords = []
    for nt in nucleotides:
        heavy_atoms = Residue3D.nucleobase_heavy_atoms.get(
            nt.one_letter_name.upper(), []
        )
        for atom_name in heavy_atoms:
            atom = nt.find_atom(atom_name)
            if atom is not None:
                coords.append(atom.coordinates)
    return coords


def calculate_rise(coords_base_A: numpy.ndarray, coords_base_B: numpy.ndarray) -> float:
    """
    Calculates the vertical rise between two lists of coordinates.

    Parameters:
    coords_base_A (numpy.ndarray): An N x 3 array of XYZ coordinates for base A.
    coords_base_B (numpy.ndarray): An M x 3 array of XYZ coordinates for base B.

    Returns:
    float: The absolute rise distance between the two bases.
    """
    # 1. Get Centroids and Normals for both bases
    c_A, n_A = get_plane_vectors(coords_base_A)
    c_B, n_B = get_plane_vectors(coords_base_B)

    # 2. Align Normals
    # If the bases are stacked, normals should point roughly the same way.
    # If dot product is negative, they point in opposite directions.
    if numpy.dot(n_A, n_B) < 0:
        n_B = -n_B

    # 3. Compute Average Normal
    # This represents the "intermediate" director axis between the two bases.
    # For the signed variant, orient the axis from A to B.
    n_avg = n_A + n_B
    n_avg = n_avg / numpy.linalg.norm(n_avg)  # Normalize to unit vector

    # 4. Compute Centroid-to-Centroid vector
    vector_connect = c_B - c_A

    # 5. Calculate Rise (Projection)
    # The dot product projects the connecting vector onto the average normal axis
    rise = numpy.dot(vector_connect, n_avg)

    # We usually care about the absolute distance, though sign indicates
    # if B is "above" or "below" A relative to the normal.
    return abs(rise)


def get_signed_angle(
    v1: numpy.ndarray, v2: numpy.ndarray, normal: numpy.ndarray
) -> float:
    """
    Calculates the signed angle (degrees) between v1 and v2,
    looking down the 'normal' axis.

    Parameters:
    v1 (numpy.ndarray): First vector.
    v2 (numpy.ndarray): Second vector.
    normal (numpy.ndarray): Normal vector defining the axis.

    Returns:
    float: Signed angle in degrees.
    """
    # 1. Project vectors onto the plane perpendicular to the normal
    # Formula: v_proj = v - (v . n) * n
    v1_proj = v1 - numpy.dot(v1, normal) * normal
    v2_proj = v2 - numpy.dot(v2, normal) * normal

    # 2. Normalize projected vectors
    v1_proj /= numpy.linalg.norm(v1_proj)
    v2_proj /= numpy.linalg.norm(v2_proj)

    # 3. Calculate Cosine (Dot Product)
    dot_prod = numpy.dot(v1_proj, v2_proj)
    # Clamp to prevent floating point errors > 1.0 or < -1.0
    dot_prod = numpy.clip(dot_prod, -1.0, 1.0)

    # 4. Calculate Sine magnitude and direction (Cross Product)
    cross_prod = numpy.cross(v1_proj, v2_proj)
    # Dotting with the normal gives the sign (Up/Down relative to normal)
    det = numpy.dot(normal, cross_prod)

    # 5. ArcTangent2 gives the full signed angle
    return numpy.degrees(numpy.arctan2(det, dot_prod))


def normalize_vector(vector: numpy.ndarray) -> numpy.ndarray:
    norm = numpy.linalg.norm(vector)
    if norm == 0.0:
        raise ValueError("Cannot normalize a zero-length vector")
    return vector / norm


def calculate_angle_around_axis(
    vector: numpy.ndarray, axis: numpy.ndarray, reference: numpy.ndarray
) -> float:
    axis = normalize_vector(axis)
    reference_proj = reference - numpy.dot(reference, axis) * axis
    vector_proj = vector - numpy.dot(vector, axis) * axis

    reference_proj = normalize_vector(reference_proj)
    vector_proj = normalize_vector(vector_proj)
    orthogonal = normalize_vector(numpy.cross(axis, reference_proj))

    return float(
        numpy.arctan2(
            numpy.dot(vector_proj, orthogonal),
            numpy.dot(vector_proj, reference_proj),
        )
    )


def calculate_best_fit_rotation_around_axis(
    points1: numpy.ndarray, points2: numpy.ndarray, axis: numpy.ndarray
) -> float:
    """Return the signed best-fit rotation between matched point sets."""
    if len(points1) != len(points2) or len(points1) == 0:
        return math.nan

    try:
        axis = normalize_vector(axis)
    except ValueError:
        return math.nan

    centered1 = points1 - numpy.mean(points1, axis=0)
    centered2 = points2 - numpy.mean(points2, axis=0)

    projected1 = centered1 - numpy.outer(centered1 @ axis, axis)
    projected2 = centered2 - numpy.outer(centered2 @ axis, axis)

    if numpy.sum(projected1 * projected1) < 1.0e-12:
        return math.nan
    if numpy.sum(projected2 * projected2) < 1.0e-12:
        return math.nan

    sine = float(numpy.sum(numpy.cross(projected1, projected2) @ axis))
    cosine = float(numpy.sum(projected1 * projected2))

    if math.isclose(sine, 0.0, abs_tol=1.0e-12) and math.isclose(
        cosine, 0.0, abs_tol=1.0e-12
    ):
        return math.nan

    return math.degrees(math.atan2(sine, cosine))


def collect_residue_base_like_atoms(residue: Residue3D) -> list[numpy.ndarray]:
    coords = []

    for atom in residue.atoms:
        if atom.name.startswith("H"):
            continue
        if atom.name in Residue3D.phosphate_atoms or atom.name in Residue3D.sugar_atoms:
            continue
        coords.append(atom.coordinates)

    return coords


def residue_topology_point(residue: Residue3D) -> numpy.ndarray | None:
    coords = collect_nucleobase_atoms((residue,))
    if not coords:
        coords = collect_residue_base_like_atoms(residue)
    if not coords:
        return None
    return numpy.mean(coords, axis=0)


def residue_base_centroid(residue: Residue3D) -> numpy.ndarray:
    point = residue_topology_point(residue)
    if point is None:
        raise ValueError(f"Missing nucleobase atoms for residue {residue.full_name}")
    return point


def residue_c1prime_coordinates(residue: Residue3D) -> numpy.ndarray | None:
    for atom_name in ("C1'", "C1*"):
        atom = residue.find_atom(atom_name)
        if atom is not None:
            return atom.coordinates
    return None


def residue_backbone_anchor(residue: Residue3D) -> numpy.ndarray:
    atom_names = ("C3'", "C4'", "C5'", "C1'", "O4'", "O5'", "O3'")
    coords = []

    for atom_name in atom_names:
        atom = residue.find_atom(atom_name)
        if atom is not None:
            coords.append(atom.coordinates)

    if coords:
        return numpy.mean(coords, axis=0)
    point = residue_topology_point(residue)
    if point is not None:
        return point
    if residue.atoms:
        return residue.atoms[0].coordinates
    raise ValueError(f"Missing anchor atoms for residue {residue.full_name}")


def residue_backbone_direction(residue: Residue3D) -> numpy.ndarray | None:
    p = residue.find_atom("P")
    c3 = residue.find_atom("C3'")
    if p is not None and c3 is not None:
        return c3.coordinates - p.coordinates
    c5 = residue.find_atom("C5'")
    if c5 is not None and c3 is not None:
        return c3.coordinates - c5.coordinates
    return None


def sugar_torsion(residue: Residue3D, atom_names: tuple[str, str, str, str]) -> float:
    atoms = [residue.find_atom(atom_name) for atom_name in atom_names]
    if any(atom is None for atom in atoms):
        return math.nan
    return torsion_angle(atoms[0], atoms[1], atoms[2], atoms[3])


def sugar_pseudorotation(residue: Residue3D) -> float:
    nu0 = sugar_torsion(residue, ("C4'", "O4'", "C1'", "C2'"))
    nu1 = sugar_torsion(residue, ("O4'", "C1'", "C2'", "C3'"))
    nu2 = sugar_torsion(residue, ("C1'", "C2'", "C3'", "C4'"))
    nu3 = sugar_torsion(residue, ("C2'", "C3'", "C4'", "O4'"))
    nu4 = sugar_torsion(residue, ("C3'", "C4'", "O4'", "C1'"))

    if any(math.isnan(nu) for nu in (nu0, nu1, nu2, nu3, nu4)):
        return math.nan

    numerator = (nu4 + nu1) - (nu3 + nu0)
    denominator = (
        2.0 * nu2 * (math.sin(math.radians(36.0)) + math.sin(math.radians(72.0)))
    )

    if math.isclose(denominator, 0.0, abs_tol=1.0e-12):
        return math.nan

    angle = math.degrees(math.atan2(numerator, denominator))
    if angle < 0.0:
        angle += 360.0
    return angle


def classify_sugar_pucker(residue: Residue3D) -> SugarPucker | None:
    phase = sugar_pseudorotation(residue)
    if math.isnan(phase):
        return None

    if phase >= 315.0 or phase < 135.0:
        return SugarPucker.NORTH
    return SugarPucker.SOUTH


def calculate_quadruplex_twist_centroids(
    tetrad1_all_coords: numpy.ndarray,
    tetrad2_all_coords: numpy.ndarray,
    nt_list1: list[numpy.ndarray],
    nt_list2: list[numpy.ndarray],
) -> float:
    """
    Calculate the legacy centroid-based twist using matched base centroids.

    Parameters:
    - tetrad1_all_coords: (N x 3) array of all atoms in tetrad 1 (for fitting the plane).
    - tetrad2_all_coords: (N x 3) array of all atoms in tetrad 2 (for fitting the plane).
    - nt_list1: List of 4 (N x 3) arrays. Each array contains coordinates for one base (a, b, c, d).
    - nt_list2: List of 4 (N x 3) arrays. Each array contains coordinates for one base (u, v, w, x).
                (Must be ordered: a stacks on u, b on v, etc.)

    Returns:
    float: The average signed twist angle in degrees.
    """
    # A. Get Global Geometries (Tetrad Centers and Common Axis)
    c_T1 = numpy.mean(tetrad1_all_coords, axis=0)
    c_T2 = numpy.mean(tetrad2_all_coords, axis=0)
    axis = calculate_tetrad_common_axis(tetrad1_all_coords, tetrad2_all_coords)
    if axis is None:
        return math.nan

    # B. Calculate Twist for each stacked pair
    twists = []

    # We assume the lists are ordered [a, b, c, d] and [u, v, w, x]
    for i in range(4):
        # Calculate Centroid of the specific nucleobase
        # (using the atoms provided in the list)
        base_centroid_1 = numpy.mean(nt_list1[i], axis=0)
        base_centroid_2 = numpy.mean(nt_list2[i], axis=0)

        # Define Radial Vectors
        # Vector from Tetrad Center -> Base Center
        v1 = base_centroid_1 - c_T1
        v2 = base_centroid_2 - c_T2

        # Calculate signed angle using the common axis
        angle = get_signed_angle(v1, v2, axis)
        twists.append(angle)

    # C. Return the average
    return numpy.mean(twists)


def calculate_tetrad_common_axis(
    tetrad1_all_coords: numpy.ndarray, tetrad2_all_coords: numpy.ndarray
) -> numpy.ndarray | None:
    c_T1, n_T1 = get_plane_vectors(tetrad1_all_coords)
    c_T2, n_T2 = get_plane_vectors(tetrad2_all_coords)

    if numpy.dot(n_T1, n_T2) < 0:
        n_T2 = -n_T2

    axis = n_T1 + n_T2
    if numpy.linalg.norm(axis) < 1.0e-6:
        return None

    axis = normalize_vector(axis)
    if numpy.dot(axis, c_T2 - c_T1) < 0:
        axis = -axis
    return axis


def calculate_quadruplex_twist_c1primes(
    tetrad1_all_coords: numpy.ndarray,
    tetrad2_all_coords: numpy.ndarray,
    tetrad1_c1primes: numpy.ndarray,
    tetrad2_c1primes: numpy.ndarray,
) -> float:
    """Calculate twist from the best-fit rotation of matched C1' points."""
    axis = calculate_tetrad_common_axis(tetrad1_all_coords, tetrad2_all_coords)
    if axis is None:
        return math.nan

    return calculate_best_fit_rotation_around_axis(
        tetrad1_c1primes,
        tetrad2_c1primes,
        axis,
    )


def generate_tetrad_label(index: int) -> str:
    if index < 0:
        raise ValueError("Tetrad index must be non-negative")

    label = ""
    index += 1

    while index > 0:
        index, remainder = divmod(index - 1, 26)
        label = chr(ord("A") + remainder) + label

    return label


@dataclass(order=True)
class Tetrad:
    @staticmethod
    def is_valid(
        nt1: Residue3D,
        nt2: Residue3D,
        nt3: Residue3D,
        nt4: Residue3D,
        pair_dictionary: dict[tuple[Residue3D, Residue3D], BasePair3D],
    ) -> bool:
        lw1 = pair_dictionary[(nt1, nt2)].lw
        lw2 = pair_dictionary[(nt2, nt3)].lw
        lw3 = pair_dictionary[(nt3, nt4)].lw
        lw4 = pair_dictionary[(nt4, nt1)].lw
        for lw_i, lw_j in ((lw1, lw4), (lw2, lw1), (lw3, lw2), (lw4, lw3)):
            if lw_i.name[1] == lw_j.name[2]:
                return False
        return True

    nt1: Residue3D
    nt2: Residue3D
    nt3: Residue3D
    nt4: Residue3D
    pair_12: BasePair3D
    pair_23: BasePair3D
    pair_34: BasePair3D
    pair_41: BasePair3D
    global_index: dict[Residue3D, int]
    onz: ONZ = field(init=False)
    gba_class: GbaTetradClassification | None = field(init=False)
    planarity_deviation: dict[str, float] = field(init=False)
    ions_channel: list[Atom] = field(default_factory=list)
    ions_outside: dict[Residue3D, list[Atom]] = field(default_factory=dict)

    def __post_init__(self):
        self.reorder_to_match_5p_3p()
        self.planarity_deviation = self.__calculate_planarity_deviation()

    def reorder_to_match_5p_3p(self):
        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = (self.global_index[nt] for nt in self.nucleotides)
        indices = sorted((ni, nj, nk, nl))
        ni, nj, nk, nl = (indices.index(x) for x in (ni, nj, nk, nl))

        nmin = min(ni, nj, nk, nl)
        if nmin == ni:
            pass
        elif nmin == nj:
            self.nt1, self.nt2, self.nt3, self.nt4 = (
                self.nt2,
                self.nt3,
                self.nt4,
                self.nt1,
            )
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_23,
                self.pair_34,
                self.pair_41,
                self.pair_12,
            )
        elif nmin == nk:
            self.nt1, self.nt2, self.nt3, self.nt4 = (
                self.nt3,
                self.nt4,
                self.nt1,
                self.nt2,
            )
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_34,
                self.pair_41,
                self.pair_12,
                self.pair_23,
            )
        else:
            self.nt1, self.nt2, self.nt3, self.nt4 = (
                self.nt4,
                self.nt1,
                self.nt2,
                self.nt3,
            )
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_41,
                self.pair_12,
                self.pair_23,
                self.pair_34,
            )

        # flip order if necessary
        if self.pair_12.score > self.pair_41.reverse.score:
            self.nt1, self.nt2, self.nt3, self.nt4 = (
                self.nt1,
                self.nt4,
                self.nt3,
                self.nt2,
            )
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_41.reverse,
                self.pair_34.reverse,
                self.pair_23.reverse,
                self.pair_12.reverse,
            )

        # ONZ and da Silva's classification are valid in 5'-3' order
        self.onz = self.__classify_onz()
        self.gba_class = self.__classify_by_gba()

    def reorder_to_match_other_tetrad(
        self, order: tuple[Residue3D, Residue3D, Residue3D, Residue3D]
    ):
        if order == (self.nt1, self.nt2, self.nt3, self.nt4):
            pass
        elif order == (self.nt2, self.nt3, self.nt4, self.nt1):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_23,
                self.pair_34,
                self.pair_41,
                self.pair_12,
            )
        elif order == (self.nt3, self.nt4, self.nt1, self.nt2):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_34,
                self.pair_41,
                self.pair_12,
                self.pair_23,
            )
        elif order == (self.nt4, self.nt1, self.nt2, self.nt3):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_41,
                self.pair_12,
                self.pair_23,
                self.pair_34,
            )
        elif order == (self.nt4, self.nt3, self.nt2, self.nt1):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_34.reverse,
                self.pair_23.reverse,
                self.pair_12.reverse,
                self.pair_41.reverse,
            )
        elif order == (self.nt3, self.nt2, self.nt1, self.nt4):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_23.reverse,
                self.pair_12.reverse,
                self.pair_41.reverse,
                self.pair_34.reverse,
            )
        elif order == (self.nt2, self.nt1, self.nt4, self.nt3):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_12.reverse,
                self.pair_41.reverse,
                self.pair_34.reverse,
                self.pair_23.reverse,
            )
        elif order == (self.nt1, self.nt4, self.nt3, self.nt2):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = (
                self.pair_41.reverse,
                self.pair_34.reverse,
                self.pair_23.reverse,
                self.pair_12.reverse,
            )
        else:
            raise RuntimeError(f"Cannot apply order: {order}")

        self.nt1, self.nt2, self.nt3, self.nt4 = order

    def __classify_onz(self) -> ONZ:
        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = (self.global_index[nt] for nt in self.nucleotides)
        indices = sorted((ni, nj, nk, nl))
        ni, nj, nk, nl = (indices.index(x) for x in (ni, nj, nk, nl))

        while ni != 0:
            ni, nj, nk, nl = nl, ni, nj, nk

        order = (nj, nk, nl)
        if order == (1, 2, 3):
            return ONZ.O_PLUS
        elif order == (3, 2, 1):
            return ONZ.O_MINUS
        elif order == (1, 3, 2):
            return ONZ.N_PLUS
        elif order == (2, 3, 1):
            return ONZ.N_MINUS
        elif order == (2, 1, 3):
            return ONZ.Z_PLUS
        elif order == (3, 1, 2):
            return ONZ.Z_MINUS

        raise RuntimeError(f"Impossible combination: {ni} {nj} {nk} {nl}")

    def __classify_by_gba(self) -> GbaTetradClassification | None:
        """
        See: Webba da Silva, M. (2007). Geometric Formalism for DNA Quadruplex Folding.
        Chemistry - A European Journal, 13(35), 9738-9745. https://doi.org/10.1002/chem.200701255

        :return: Classification according to Webba da Silva or n/a
        """
        # without all nucleotides having a valid syn/anti, this classification is impossible
        if not all(
            nt.chi_class in (GlycosidicBond.syn, GlycosidicBond.anti)
            for nt in self.nucleotides
        ):
            return None

        # this will create a 4-letter string made of 's' for syn or 'a' for anti
        fingerprint = "".join([nt.chi_class.value[0] for nt in self.nucleotides])  # type: ignore

        # this dict has all classes mapped to fingerprints
        gba_classes = {
            "aass": GbaTetradClassification.Ia,
            "ssaa": GbaTetradClassification.Ib,
            "asas": GbaTetradClassification.IIa,
            "sasa": GbaTetradClassification.IIb,
            "asaa": GbaTetradClassification.IIIa,
            "sass": GbaTetradClassification.IIIb,
            "aaas": GbaTetradClassification.IVa,
            "sssa": GbaTetradClassification.IVb,
            "aasa": GbaTetradClassification.Va,
            "ssas": GbaTetradClassification.Vb,
            "assa": GbaTetradClassification.VIa,
            "saas": GbaTetradClassification.VIb,
            "asss": GbaTetradClassification.VIIa,
            "saaa": GbaTetradClassification.VIIb,
            "aaaa": GbaTetradClassification.VIIIa,
            "ssss": GbaTetradClassification.VIIIb,
        }

        if fingerprint not in gba_classes:
            logger.error(
                "Impossible combination of syn/anti: %s",
                [nt.chi_class for nt in self.nucleotides],
            )
            return None
        return gba_classes[fingerprint]

    def __calculate_planarity_deviation(self) -> dict[str, float]:
        atom_coords = collect_nucleobase_atoms(self.nucleotides)
        return calculate_planarity_deviation(numpy.array(atom_coords))

    @property
    def nucleotides(self) -> tuple[Residue3D, Residue3D, Residue3D, Residue3D]:
        return self.nt1, self.nt2, self.nt3, self.nt4

    def __hash__(self):
        return hash(frozenset([self.nt1, self.nt2, self.nt3, self.nt4]))

    def __repr__(self):
        return f"{self.nt1!r}-{self.nt2!r}-{self.nt3!r}-{self.nt4!r}"

    def __str__(self):
        return (
            f"    "
            f"{self.nt1.full_name} {self.nt2.full_name} {self.nt3.full_name} {self.nt4.full_name} "
            f"{self.pair_12.lw.value} {self.pair_23.lw.value} {self.pair_34.lw.value} {self.pair_41.lw.value} "
            f"{self.onz.value} {self.gba_class.value if self.gba_class is not None else ''} "
            f"planarity={round(self.planarity_deviation['rmsd'], 2)} "
            f"{self.__ions_channel_str()} "
            f"{self.__ions_outside_str()}\n"
        )

    def chains(self) -> set[str]:
        return {nt.chain for nt in self.nucleotides}

    def is_disjoint(self, other) -> bool:
        return frozenset(self.nucleotides).isdisjoint(frozenset(other.nucleotides))

    def center(self) -> numpy.ndarray:
        return center_of_mass(self.outer_and_inner_atoms())

    def outer_and_inner_atoms(self) -> list[Atom]:
        return [residue.outermost_atom for residue in self.nucleotides] + [
            residue.innermost_atom for residue in self.nucleotides
        ]

    def __ions_channel_str(self) -> str:
        if self.ions_channel:
            return "ions_channel=" + ",".join([atom.name for atom in self.ions_channel])
        return ""

    def __ions_outside_str(self) -> str:
        if self.ions_outside:
            result = []
            for residue, ions in self.ions_outside.items():
                result.append(
                    f"{residue.full_name}: [{','.join([ion.name for ion in ions])}]"
                )
            return "ions_outside=" + " ".join(result)
        return ""


@dataclass
class TetradPair:
    """Describe one adjacent stacked tetrad pair.

    ``tetrad2_nts_best_order`` reorders the second tetrad so residue ``i`` in the
    first tetrad is the residue that best continues the same tract in the second
    tetrad. The reported ``rise`` and ``twist`` therefore describe adjacent stack
    geometry, not arbitrary build-order intervals used by g4composer.
    """

    tetrad1: Tetrad
    tetrad2: Tetrad
    stacked: dict[Residue3D, Residue3D]
    global_index: dict[Residue3D, int]
    tetrad2_nts_best_order: tuple[Residue3D, Residue3D, Residue3D, Residue3D] = field(
        init=False
    )
    direction: Direction = field(init=False)
    rise: float = field(init=False)
    twist: float = field(init=False)

    def __post_init__(self):
        self.tetrad2_nts_best_order = (
            self.stacked[self.tetrad1.nt1],
            self.stacked[self.tetrad1.nt2],
            self.stacked[self.tetrad1.nt3],
            self.stacked[self.tetrad1.nt4],
        )
        self.direction = self.__determine_direction()
        self.rise = self.__calculate_rise()
        self.twist = self.__calculate_twist()

    def __determine_direction(self) -> Direction:
        indices1 = [self.global_index[nt] for nt in self.tetrad1.nucleotides]
        indices2 = [self.global_index[nt] for nt in self.tetrad2_nts_best_order]

        # count directions 5' -> 3' as +1 or -1
        counter = Counter(
            1 if j - i > 0 else -1 for i, j in zip(indices1, indices2, strict=True)
        )
        _direction, count = counter.most_common()[0]

        if count == 4:
            # all in the same direction
            return Direction.parallel
        elif count == 2:
            # two in +, one in - direction
            return Direction.antiparallel

        return Direction.hybrid

    def __calculate_rise(self) -> float:
        coords_base_A = collect_nucleobase_atoms(self.tetrad1.nucleotides)
        coords_base_B = collect_nucleobase_atoms(self.tetrad2.nucleotides)

        if coords_base_A and coords_base_B:
            return calculate_rise(
                numpy.array(coords_base_A), numpy.array(coords_base_B)
            )
        return math.nan

    def __calculate_twist(self) -> float:
        """Measure twist from the best-fit rotation of matched C1' atoms."""
        # Collect all nucleobase atoms for tetrad 1
        tetrad1_all_coords = collect_nucleobase_atoms(self.tetrad1.nucleotides)
        # Collect all nucleobase atoms for tetrad 2
        tetrad2_all_coords = collect_nucleobase_atoms(self.tetrad2_nts_best_order)

        if not tetrad1_all_coords or not tetrad2_all_coords:
            return math.nan

        tetrad1_c1primes = []
        for nt in self.tetrad1.nucleotides:
            c1prime = residue_c1prime_coordinates(nt)
            if c1prime is None:
                return math.nan
            tetrad1_c1primes.append(c1prime)

        tetrad2_c1primes = []
        for nt in self.tetrad2_nts_best_order:
            c1prime = residue_c1prime_coordinates(nt)
            if c1prime is None:
                return math.nan
            tetrad2_c1primes.append(c1prime)

        return calculate_quadruplex_twist_c1primes(
            numpy.array(tetrad1_all_coords),
            numpy.array(tetrad2_all_coords),
            numpy.array(tetrad1_c1primes),
            numpy.array(tetrad2_c1primes),
        )

    def __str__(self):
        return f"      direction={self.direction.value} rise={round(self.rise, 2)} twist={round(self.twist, 2)}\n"


@dataclass
@dataclass(frozen=True)
class TetradScore:
    """
    Detailed description of the alignment between two tetrads.
    """

    sequential: int  # strictly sequential nucleotide matches
    stacking: int  # stacking-graph matches
    nts1: tuple[Residue3D, ...]  # order in the first tetrad
    nts2: tuple[Residue3D, ...]  # corresponding order in the second tetrad


@dataclass
class Tract:
    nucleotides: list[Residue3D]

    def __str__(self):
        return f"      {', '.join(nt.full_name for nt in self.nucleotides)}"


@dataclass
class Loop:
    nucleotides: list[Residue3D]
    loop_type: LoopType | None
    missing: list[str] = field(default_factory=list)

    def __str__(self):
        return (
            f"      {self.loop_type.value if self.loop_type else 'n/a'} "
            f"{', '.join(nt.full_name for nt in self.nucleotides)}"
            f"{self.__missing_suffix()}"
        )

    def __missing_suffix(self) -> str:
        return f" (missing: {', '.join(self.missing)})" if self.missing else ""


@dataclass
class Snapback:
    """A chain re-entry into an already visited tract (same column, tetrad
    layers not adjacent). The nucleotides are the connector between the two
    same-column Gs (possibly empty); orientation says which side of the
    connector is the isolated one."""

    nucleotides: list[Residue3D]
    orientation: SnapbackOrientation
    missing: list[str] = field(default_factory=list)

    def __str__(self):
        return (
            f"      snapback {self.orientation.value} "
            f"{', '.join(nt.full_name for nt in self.nucleotides)}"
            f"{self.__missing_suffix()}"
        )

    def __missing_suffix(self) -> str:
        return f" (missing: {', '.join(self.missing)})" if self.missing else ""


@dataclass
class BulgeGroup:
    """Nucleotides bulging out between two stacked nucleotides of one tract,
    together with any residues reported missing between the flanks."""

    nucleotides: list[Residue3D]
    missing: list[str] = field(default_factory=list)

    def __str__(self):
        names = [nt.full_name for nt in self.nucleotides]
        names.extend(self.missing)
        return ", ".join(names)


# Mapping from the Webba da Silva loop-classification value (the loop
# progression fingerprint with the a/b subtype given by the sign of the
# first flanking loop, or of the middle loop for the dpd progression) to
# the common topology name.
#
# Source: ASC-G4, "Comprehensive analysis of intramolecular G-quadruplex
# structures: furthering the understanding of their formalism"
# (Nucleic Acids Res. 2024, gkae182; PMC11039995), Table 1, which lists the
# observed loop progressions per topology in the Dvorkin et al. (2018)
# frame of reference (the same column anchoring ElTetrado uses):
#   parallel  -p-p-p                 (1a; 1b +p+p+p observed in PDB as well)
#   chair     +l+l+l, -l-l-l         (6b, 6a)
#   basket    -ld+l, d+pd            (11a, 12b)
#   basket2   -pd+p, +ld-l, +l+p+l   (5a, 11b, 8b)
#   hybrid1   -p-l-l, -p-l-p         (9a, 3a)
#   hybrid2   -pd+l, -p-p-l, +ld-p   (10a, 2a, 13b)
#   hybrid3   -l-l-p                 (7a)
#   hybrid4   +l+p+p                 (4b)
# The remaining subtypes (12a, 5b, 8a, 3b, 9b, 2b, 10b, 13a, 7b, 4a) have
# not been observed and carry no common name - see AGENTS.md.
TOPOLOGY_NAMES = {
    "1a": "parallel",
    "1b": "parallel",
    "6a": "chair",
    "6b": "chair",
    "11a": "basket",
    "12b": "basket",
    "5a": "basket2",
    "11b": "basket2",
    "8b": "basket2",
    "9a": "hybrid1",
    "3a": "hybrid1",
    "10a": "hybrid2",
    "2a": "hybrid2",
    "13b": "hybrid2",
    "7a": "hybrid3",
    "4b": "hybrid4",
}


def sequence_entries_with_missing(
    nucleotides: list[Residue3D],
    missing_residues: "list[MissingResidue] | None",
    sequence_by_entity: "dict[str, str] | None",
) -> list[tuple[str, str, Residue3D | None]]:
    """Build the full per-chain sequence entries for the report and exports.

    Returns (chain, one-letter code, residue or None) triples in the given
    nucleotide order (already chain-ordered by the caller), with residues
    reported missing in the file header interleaved at their numbering
    positions. Missing residues take their real one-letter code from the
    entity sequence when available and carry ``residue=None``; without the
    header data only numbering-based '?' placeholders (unknown identity) are
    inserted for internal breaks.
    """
    if not missing_residues:
        # degraded fallback: only internal numbering gaps between
        # geometrically disconnected residues are marked with '?'
        entries: list[tuple[str, str, Residue3D | None]] = []
        for position, residue in enumerate(nucleotides):
            if (
                position > 0
                and residue.chain == nucleotides[position - 1].chain
                and not nucleotides[position - 1].is_connected(residue)
            ):
                previous_number = nucleotides[position - 1].number
                modelled_numbers = {
                    other.number
                    for other in nucleotides
                    if other.chain == residue.chain
                }
                for number in range(
                    min(previous_number, residue.number) + 1,
                    max(previous_number, residue.number),
                ):
                    if number not in modelled_numbers:
                        entries.append((residue.chain, "?", None))
            entries.append((residue.chain, residue.one_letter_name, residue))
        return entries

    modelled_by_chain: dict[str, dict[int, Residue3D]] = {}
    chain_order: list[str] = []
    entity_by_chain: dict[str, str | None] = {}
    for residue in nucleotides:
        modelled_by_chain.setdefault(residue.chain, {})[residue.number] = residue
        if residue.chain not in chain_order:
            chain_order.append(residue.chain)
            entity_by_chain[residue.chain] = (
                residue.atoms[0].entity_id if residue.atoms else None
            )

    missing_by_chain: dict[str, dict[int, MissingResidue]] = {}
    for missing in missing_residues:
        if missing.chain not in modelled_by_chain or missing.number is None:
            continue
        missing_by_chain.setdefault(missing.chain, {})[missing.number] = missing

    entries: list[tuple[str, str, Residue3D | None]] = []
    for chain in chain_order:
        modelled = modelled_by_chain[chain]
        missing = missing_by_chain.get(chain, {})
        entity_id = entity_by_chain.get(chain)
        sequence = sequence_by_entity.get(entity_id) if entity_id else None
        for number in sorted(set(modelled) | set(missing)):
            residue = modelled.get(number)
            if residue is not None:
                entries.append((chain, residue.one_letter_name, residue))
                continue
            missing_residue = missing[number]
            letter = missing_residue.one_letter_name
            label_number = (
                missing_residue.label.number if missing_residue.label else None
            )
            if (
                sequence is not None
                and label_number
                and 1 <= label_number <= len(sequence)
            ):
                letter = sequence[label_number - 1]
            entries.append((chain, letter, None))
    return entries


def quadruplex_spurious_reasons(quadruplex: "Quadruplex") -> list[str]:
    """Return reasons why a quadruplex should be discarded as spurious.

    An empty list means the quadruplex passed all filter limits.
    """
    reasons = []

    if len(quadruplex.tetrads) < MIN_QUADRUPLEX_TETRADS:
        reasons.append(
            f"built of {len(quadruplex.tetrads)} tetrad(s), "
            f"minimum is {MIN_QUADRUPLEX_TETRADS}"
        )

    bulge_count = sum(len(bulge.nucleotides) for bulge in quadruplex.bulges)
    if bulge_count > MAX_QUADRUPLEX_BULGES:
        reasons.append(
            f"bulge of {bulge_count} nucleotides (limit is {MAX_QUADRUPLEX_BULGES})"
        )

    for tetrad_pair in quadruplex.tetrad_pairs:
        rise = tetrad_pair.rise
        if not math.isnan(rise) and not (
            MIN_QUADRUPLEX_RISE <= rise <= MAX_QUADRUPLEX_RISE
        ):
            reasons.append(
                f"rise between tetrads {round(rise, 2)} Å "
                f"(allowed range is {MIN_QUADRUPLEX_RISE}-{MAX_QUADRUPLEX_RISE} Å)"
            )

    return reasons


def is_pure_g_tetrad(tetrad: "Tetrad") -> bool:
    """True when all tetrad nucleotides are guanines (including modified
    guanines such as 8OG or LCG)."""
    return all(nt.one_letter_name.upper() == "G" for nt in tetrad.nucleotides)


def quadruplex_tract_runs(quadruplex: "Quadruplex") -> list[list[list[Residue3D]]]:
    """Split each tract into runs of nucleotides that are consecutive in the
    chain-ordered sequence of all tetrad nucleotides (result is aligned with
    ``quadruplex.tracts``).

    A run breaks wherever the chain leaves the tract and visits another
    tract's tetrad nucleotide(s) before coming back (a re-entry excursion).
    Bulges between two stacked nucleotides of the same tract and direct 0-nt
    reversals do NOT break a run.
    """
    tract_of: dict[Residue3D, int] = {}
    for tract_index, tract in enumerate(quadruplex.tracts):
        for nt in tract.nucleotides:
            if nt not in tract_of:
                tract_of[nt] = tract_index
    chain_ordered = sorted(tract_of, key=lambda nt: quadruplex.global_index[nt])
    runs: list[list[list[Residue3D]]] = [[] for _ in quadruplex.tracts]
    last_seen_at: dict[int, int] = {}
    last_run: dict[int, list[Residue3D]] = {}
    for position, nt in enumerate(chain_ordered):
        tract_index = tract_of[nt]
        if (
            last_seen_at.get(tract_index) == position - 1
            and last_run[tract_index][-1].chain == nt.chain
        ):
            last_run[tract_index].append(nt)
        else:
            last_run[tract_index] = [nt]
            runs[tract_index].append(last_run[tract_index])
        last_seen_at[tract_index] = position
    return runs


def isolated_run_nucleotides(
    tract_runs: list[list[list[Residue3D]]],
) -> set[Residue3D]:
    """Nucleotides forming a single-nucleotide run of a tract that the chain
    visits in more than one run (the "snapped" guanines of re-entries)."""
    isolated: set[Residue3D] = set()
    for runs in tract_runs:
        if len(runs) >= 2:
            for run in runs:
                if len(run) == 1:
                    isolated.add(run[0])
    return isolated


def quadruplex_row_partition_boundaries(quadruplex: "Quadruplex") -> list[int]:
    """Tetrad-row boundaries at which every column splits into two
    chain-homogeneous halves (each side of the boundary is owned by a single
    chain in every tract). Empty for staggered arrangements such as 1jjp or
    1nyd; vacuously full for unimolecular quadruplexes."""
    tetrad_count = len(quadruplex.tetrads)
    return [
        boundary
        for boundary in range(1, tetrad_count)
        if all(
            len({nt.chain for nt in tract.nucleotides[:boundary]}) <= 1
            and len({nt.chain for nt in tract.nucleotides[boundary:]}) <= 1
            for tract in quadruplex.tracts
        )
    ]


@dataclass
class Quadruplex:
    tetrads: list[Tetrad]
    tetrad_pairs: list[TetradPair]
    structure3d: Structure3D
    global_index: dict[Residue3D, int]
    base_pairs: list[BasePair3D] = field(default_factory=list)
    onzm: ONZM | None = field(init=False)
    gba_classes: list[GbaQuadruplexClassification] = field(init=False)
    tracts: list[Tract] = field(init=False)
    bulges: list[Residue3D] = field(init=False)
    loops: list[Loop] = field(init=False)
    snapbacks: list[Snapback] = field(init=False)
    tags: list[str] = field(init=False)
    path: list[str] = field(init=False)
    handedness: HelixHandedness | None = field(init=False)
    tetrad_polarities: list[TetradPolarity | None] = field(init=False)
    strand_polarities: list[list[StrandPolarity | None]] = field(init=False)
    loop_class: LoopClassification | None = field(init=False)

    def __post_init__(self):
        self.onzm = self.__classify_onzm()
        self.gba_classes = self.__classify_by_gba()
        self.tracts = self.__find_tracts()
        self.bulges = self.__find_bulges()
        self.loops = self.__find_loops()
        self.snapbacks = self.__find_snapbacks()
        self.handedness = self.__find_handedness()
        self.tetrad_polarities = self.__find_tetrad_polarities()
        self.strand_polarities = self.__find_strand_polarities()
        self.path = self.__find_path()
        self.loop_class = self.__classify_by_loops()
        self.tags = self.__find_tags()

    def __classify_onzm(self) -> ONZM | None:
        if len(self.tetrads) == 1:
            return None
        if any(t.onz is None for t in self.tetrads):
            return None

        counter = Counter([t.onz.value[0] for t in self.tetrads])
        onz, support = counter.most_common()[0]
        if support != len(self.tetrads):
            onz = "M"

        counter = Counter([tp.direction.value[0] for tp in self.tetrad_pairs])
        direction, support = counter.most_common()[0]
        if support != len(self.tetrad_pairs):
            direction = "h"

        counter = Counter([t.onz.value[1] for t in self.tetrads])
        plus_minus, support = counter.most_common()[0]
        if support != len(self.tetrads):
            plus_minus = "*"

        return ONZM.from_value(f"{onz}{direction}{plus_minus}")

    def __classify_by_gba(self) -> list[GbaQuadruplexClassification]:
        gbas = set()
        for t in self.tetrads:
            gba = t.gba_class
            if gba is not None:
                gbas.add(gba.value[:-1])  # discard 'a' or 'b' subvariant
        roman_numerals = {
            "I": 1,
            "II": 2,
            "III": 3,
            "IV": 4,
            "V": 5,
            "VI": 6,
            "VII": 7,
            "VIII": 8,
        }
        return [
            GbaQuadruplexClassification[x]
            for x in sorted(gbas, key=lambda gba: roman_numerals.get(gba, 100))
        ]

    def __find_tracts(self) -> list[Tract]:
        tracts = [
            Tract(nts)
            for nts in self.__build_tracts(self.__first_tetrad_column_order())
        ]
        return self.__rotate_tracts_to_tetrad_a(tracts)

    def __rotate_tracts_to_tetrad_a(self, tracts: list[Tract]) -> list[Tract]:
        """Re-anchor tract numbering so tract ``1`` always belongs to tetrad A.

        ``__first_tetrad_column_order`` anchors column ``1`` on
        ``self.tetrads[0]``, which is the first tetrad in the internal 3D-stack
        discovery order, not necessarily the tetrad labelled ``A`` (the tetrad
        with the 5'-most nucleotide). Rotating the tract list preserves the
        existing geometric (anticlockwise) cyclic order while guaranteeing that
        ``path`` always starts with ``A1``.
        """
        if not tracts or not self.tetrads:
            return tracts

        tetrad_a_index = min(
            range(len(self.tetrads)),
            key=lambda i: min(
                self.global_index[nt] for nt in self.tetrads[i].nucleotides
            ),
        )
        tetrad_a_anchor = min(
            self.tetrads[tetrad_a_index].nucleotides,
            key=lambda nt: self.global_index[nt],
        )

        for rotation, tract in enumerate(tracts):
            if tract.nucleotides[tetrad_a_index] == tetrad_a_anchor:
                return tracts[rotation:] + tracts[:rotation]

        return tracts

    def __build_tracts(
        self, first_tetrad_order: list[Residue3D]
    ) -> list[list[Residue3D]]:
        tracts = [[nt] for nt in first_tetrad_order]

        for tetrad_pair in self.tetrad_pairs:
            nt_dict = {
                tetrad_pair.tetrad1.nt1: tetrad_pair.tetrad2_nts_best_order[0],
                tetrad_pair.tetrad1.nt2: tetrad_pair.tetrad2_nts_best_order[1],
                tetrad_pair.tetrad1.nt3: tetrad_pair.tetrad2_nts_best_order[2],
                tetrad_pair.tetrad1.nt4: tetrad_pair.tetrad2_nts_best_order[3],
            }
            for tract in tracts:
                tract.append(nt_dict[tract[-1]])

        return tracts

    def __first_tetrad_column_order(self) -> list[Residue3D]:
        tetrad = self.tetrads[0]
        geometry = self.__tetrad_geometry(tetrad)
        if geometry is None:
            return list(tetrad.nucleotides)

        center, normal = geometry
        first_nt = min(tetrad.nucleotides, key=lambda nt: self.global_index[nt])
        normal = self.__orient_tetrad_normal(0, tetrad, normal, first_nt)

        reference_point = residue_topology_point(first_nt)
        if reference_point is None:
            return list(tetrad.nucleotides)

        reference = reference_point - center
        angles = {}
        try:
            for nt in tetrad.nucleotides:
                point = residue_topology_point(nt)
                if point is None:
                    return list(tetrad.nucleotides)
                angles[nt] = calculate_angle_around_axis(
                    point - center,
                    normal,
                    reference,
                )
        except ValueError:
            return list(tetrad.nucleotides)

        anchor_angle = angles[first_nt]

        return sorted(
            tetrad.nucleotides,
            key=lambda nt: (angles[nt] - anchor_angle) % (2.0 * math.pi),
        )

    def __tetrad_geometry(
        self, tetrad: Tetrad
    ) -> tuple[numpy.ndarray, numpy.ndarray] | None:
        atom_coords = collect_nucleobase_atoms(tetrad.nucleotides)
        if len(atom_coords) >= 3:
            return get_plane_vectors(numpy.array(atom_coords))

        points = [
            point
            for point in (residue_topology_point(nt) for nt in tetrad.nucleotides)
            if point is not None
        ]
        if len(points) >= 3:
            return get_plane_vectors(numpy.array(points))

        return None

    def __orient_tetrad_normal(
        self,
        tetrad_index: int,
        tetrad: Tetrad,
        normal: numpy.ndarray,
        reference_nt: Residue3D,
    ) -> numpy.ndarray:
        direction = self.__strand_progression_vector(tetrad_index, tetrad, reference_nt)
        if direction is not None and numpy.dot(normal, direction) < 0:
            return -normal
        return normal

    def __strand_progression_vector(
        self, tetrad_index: int, tetrad: Tetrad, residue: Residue3D
    ) -> numpy.ndarray | None:
        if tetrad_index < len(self.tetrad_pairs):
            stacked = self.tetrad_pairs[tetrad_index].stacked.get(residue)
            if stacked is not None:
                return residue_backbone_anchor(stacked) - residue_backbone_anchor(
                    residue
                )

        if tetrad_index > 0:
            stacked = self.tetrad_pairs[tetrad_index - 1].stacked.get(residue)
            if stacked is not None:
                return residue_backbone_anchor(residue) - residue_backbone_anchor(
                    stacked
                )

        return self.__local_backbone_direction(residue)

    def __local_backbone_direction(self, residue: Residue3D) -> numpy.ndarray | None:
        nucleotides = [
            nt
            for nt in self.structure3d.residues
            if nt.is_nucleotide and nt.chain == residue.chain
        ]
        nucleotides.sort(key=lambda nt: self.global_index[nt])

        try:
            position = nucleotides.index(residue)
        except ValueError:
            return None

        if position + 1 < len(nucleotides):
            return residue_backbone_anchor(
                nucleotides[position + 1]
            ) - residue_backbone_anchor(residue)
        if position > 0:
            return residue_backbone_anchor(residue) - residue_backbone_anchor(
                nucleotides[position - 1]
            )

        atom_o3 = residue.find_atom("O3'")
        atom_o5 = residue.find_atom("O5'")
        if atom_o3 is not None and atom_o5 is not None:
            return atom_o3.coordinates - atom_o5.coordinates

        return None

    def __chain_connectors(self):
        """Yield ``(nt_prev, nt_cur, in_between)`` for every pair of
        chain-consecutive tetrad nucleotides (same chain, in chain order).
        ``in_between`` holds the standard nucleotides strictly between the
        two and may be empty. This single walk classifies every connector of
        the chain through the tetrad core: cross-tract connectors are loops,
        same-tract connectors between adjacent tract positions are stacks
        (bulging out any in-between nucleotides) and same-tract connectors
        between non-adjacent tract positions are snapbacks."""
        tetrad_nucleotides = sorted(
            [nt for tetrad in self.tetrads for nt in tetrad.nucleotides],
            key=lambda nt: self.global_index[nt],
        )
        for i in range(1, len(tetrad_nucleotides)):
            nt_prev = tetrad_nucleotides[i - 1]
            nt_cur = tetrad_nucleotides[i]
            if nt_prev.chain != nt_cur.chain:
                continue
            in_between = list(
                filter(
                    lambda nt: (
                        nt.is_nucleotide
                        and self.global_index[nt_prev]
                        < self.global_index[nt]
                        < self.global_index[nt_cur]
                    ),
                    self.structure3d.residues,
                )
            )
            yield nt_prev, nt_cur, in_between

    def __find_loops(self) -> list[Loop]:
        if len(self.tetrads) == 1:
            return []

        loops = []
        for nt_prev, nt_cur, in_between in self.__chain_connectors():
            tract_prev = self.__find_tract_index_with_nt(nt_prev)
            tract_cur = self.__find_tract_index_with_nt(nt_cur)
            if tract_prev is None or tract_cur is None or tract_prev == tract_cur:
                continue
            loop_type = self.__detect_loop_type(nt_prev, nt_cur)
            loops.append(
                Loop(
                    in_between,
                    loop_type,
                    self.__missing_residue_names_between(nt_prev, nt_cur),
                )
            )
        return loops

    def __find_path(self) -> list[str]:
        labels = {}
        tetrad_labels = self.__tetrad_labels_by_5p_order()

        for tetrad_index, tetrad_label in enumerate(tetrad_labels):
            for tract_index, tract in enumerate(self.tracts, start=1):
                labels[tract.nucleotides[tetrad_index]] = f"{tetrad_label}{tract_index}"

        return [
            labels[nt]
            for nt in sorted(labels.keys(), key=lambda nt: self.global_index[nt])
        ]

    def __tetrad_labels_by_5p_order(self) -> list[str]:
        tetrad_indices = sorted(
            range(len(self.tetrads)),
            key=lambda i: min(
                self.global_index[nt] for nt in self.tetrads[i].nucleotides
            ),
        )
        labels = [""] * len(self.tetrads)

        for label_index, tetrad_index in enumerate(tetrad_indices):
            labels[tetrad_index] = generate_tetrad_label(label_index)

        return labels

    def tetrad_labels_by_5p_order(self) -> list[str]:
        return self.__tetrad_labels_by_5p_order()

    def tetrad_geometry(
        self, tetrad: Tetrad
    ) -> tuple[numpy.ndarray, numpy.ndarray] | None:
        return self.__tetrad_geometry(tetrad)

    def oriented_tetrad_normal(
        self,
        tetrad_index: int,
        tetrad: Tetrad,
        normal: numpy.ndarray,
        reference_nt: Residue3D,
    ) -> numpy.ndarray:
        return self.__orient_tetrad_normal(tetrad_index, tetrad, normal, reference_nt)

    def __find_handedness(self) -> HelixHandedness | None:
        if len(self.tetrads) < 2:
            return None

        axis = self.__compute_axis()
        if axis is None:
            return None

        centers = [tetrad.center() for tetrad in self.tetrads]

        twists = []
        for i in range(1, len(self.tetrads)):
            center_prev = centers[i - 1]
            center_cur = centers[i]
            for tract in self.tracts:
                prev_anchor = residue_backbone_anchor(tract.nucleotides[i - 1])
                cur_anchor = residue_backbone_anchor(tract.nucleotides[i])
                twists.append(
                    get_signed_angle(
                        prev_anchor - center_prev,
                        cur_anchor - center_cur,
                        axis,
                    )
                )

        if not twists:
            return None

        mean_twist = float(numpy.mean(twists))
        if math.isclose(mean_twist, 0.0, abs_tol=1.0e-6):
            return None
        return HelixHandedness.RIGHT if mean_twist > 0.0 else HelixHandedness.LEFT

    def __compute_axis(self) -> numpy.ndarray | None:
        normals = []
        for tetrad in self.tetrads:
            geometry = self.__tetrad_geometry(tetrad)
            if geometry is None:
                continue
            _center, normal = geometry
            normal = self.__orient_tetrad_normal(
                self.tetrads.index(tetrad),
                tetrad,
                normal,
                self.tracts[0].nucleotides[0],
            )
            normals.append(normal)

        if not normals:
            return None

        # Ensure consistent orientation of all normals
        for i in range(1, len(normals)):
            if numpy.dot(normals[0], normals[i]) < 0:
                normals[i] = -normals[i]

        axis = numpy.mean(normals, axis=0)
        if numpy.linalg.norm(axis) == 0.0:
            return None
        axis = normalize_vector(axis)

        # Orient the axis using the first tract's backbone progression
        if self.tracts:
            progression = self.__strand_progression_vector(
                0, self.tetrads[0], self.tracts[0].nucleotides[0]
            )
            if progression is not None and numpy.dot(axis, progression) < 0:
                axis = -axis
        return axis

    def __find_strand_polarities(self) -> list[list[StrandPolarity | None]]:
        axis = self.__compute_axis()
        if axis is None:
            return [[None] * len(t.nucleotides) for t in self.tracts]

        # Compute dot products for all nucleotides in all tetrads
        all_dots: dict[Residue3D, float | None] = {}
        for tract in self.tracts:
            for nt in tract.nucleotides:
                v = residue_backbone_direction(nt)
                if v is not None:
                    all_dots[nt] = float(numpy.dot(v, axis))
                else:
                    all_dots[nt] = None

        # Find the most 5' guanine with a clear dot product as the reference
        threshold = 0.1
        reference_dot: float | None = None
        for nt in sorted(all_dots.keys(), key=lambda n: self.global_index[n]):
            d = all_dots[nt]
            if d is not None and abs(d) > threshold:
                reference_dot = d
                break

        if reference_dot is None:
            return [[None] * len(t.nucleotides) for t in self.tracts]

        result = []
        for tract in self.tracts:
            tract_polarities = []
            for nt in tract.nucleotides:
                d = all_dots.get(nt)
                if d is None:
                    tract_polarities.append(None)
                elif d * reference_dot > threshold * threshold:
                    # Same direction as reference = normal/upright
                    tract_polarities.append(StrandPolarity.PLUS)
                elif d * reference_dot < -threshold * threshold:
                    # Opposite direction to reference = inverted
                    tract_polarities.append(StrandPolarity.MINUS)
                else:
                    tract_polarities.append(None)
            result.append(tract_polarities)
        return result

    def __find_tetrad_polarities(self) -> list[TetradPolarity | None]:
        """Judge each tetrad's H-bond successor cycle around the g4 axis.

        The g4 axis is the analysis axis re-oriented so that the centroid
        displacement between the first two tetrads in build (path) order is
        positive - the direction g4composer describes the quadruplex in. The
        polarity is CLOCKWISE when the hydrogen-bond successor cycle runs in
        the direction that is clockwise when viewed along that axis
        (guanines ordered by decreasing angle around it), ANTICLOCKWISE when
        it runs the opposite way.
        """
        axis = self.__polarities_axis()
        polarities: list[TetradPolarity | None] = []

        for tetrad in self.tetrads:
            geometry = self.__tetrad_geometry(tetrad)
            if geometry is None or axis is None:
                polarities.append(None)
                continue

            center, _ = geometry
            successor = self.__hydrogen_bond_successor_map(tetrad)
            if len(successor) != 4:
                polarities.append(None)
                continue

            ordered = self.__axis_clockwise_order(tetrad, center, axis)
            if ordered is None:
                polarities.append(None)
                continue

            next_nt = successor.get(ordered[0])
            if next_nt == ordered[1]:
                polarities.append(TetradPolarity.CLOCKWISE)
            elif next_nt == ordered[3]:
                polarities.append(TetradPolarity.ANTICLOCKWISE)
            else:
                polarities.append(None)

        return polarities

    def __polarities_axis(self) -> numpy.ndarray | None:
        """The g4 axis for the polarity judgement: the analysis axis oriented
        so that the first build-order interval's centroid displacement is
        positive."""
        axis = self.__compute_axis()
        if axis is None:
            return None
        build_order = self.__build_order_tetrads()
        if len(build_order) < 2:
            return axis
        first = self.__tetrad_geometry(build_order[0])
        second = self.__tetrad_geometry(build_order[1])
        if first is None or second is None:
            return axis
        displacement = second[0] - first[0]
        if numpy.dot(axis, displacement) < 0:
            axis = -axis
        return axis

    def __build_order_tetrads(self) -> list[Tetrad]:
        """Tetrads in g4composer build order: the first-appearance order of
        the tetrad labels along the chain."""
        labels = self.tetrad_labels_by_5p_order()
        label_of: dict[Residue3D, str] = {
            nt: labels[index]
            for index, tetrad in enumerate(self.tetrads)
            for nt in tetrad.nucleotides
        }
        ordered: list[Tetrad] = []
        seen: set[str] = set()
        for nt in sorted(label_of, key=lambda nt: self.global_index[nt]):
            label = label_of[nt]
            if label in seen:
                continue
            seen.add(label)
            ordered.append(
                next(tetrad for tetrad in self.tetrads if nt in tetrad.nucleotides)
            )
        return ordered

    def __axis_clockwise_order(
        self,
        tetrad: Tetrad,
        center: numpy.ndarray,
        axis: numpy.ndarray,
    ) -> list[Residue3D] | None:
        """The tetrad's nucleotides ordered clockwise when viewed along the
        g4 axis: by decreasing angle of C1' (or the topology point) projected
        onto the plane perpendicular to the axis."""
        angles: list[tuple[float, Residue3D]] = []
        for nt in tetrad.nucleotides:
            point = residue_c1prime_coordinates(nt)
            if point is None:
                point = residue_topology_point(nt)
            if point is None:
                return None
            vector = point - center
            vector = vector - numpy.dot(vector, axis) * axis
            reference = numpy.array([1.0, 0.0, 0.0])
            if abs(float(numpy.dot(reference, axis))) > 0.9:
                reference = numpy.array([0.0, 1.0, 0.0])
            reference = reference - numpy.dot(reference, axis) * axis
            angle = math.atan2(
                float(numpy.dot(numpy.cross(reference, vector), axis)),
                float(numpy.dot(reference, vector)),
            )
            angles.append((angle, nt))
        angles.sort(key=lambda item: item[0], reverse=True)
        return [nt for _, nt in angles]

    def __hydrogen_bond_successor_map(
        self, tetrad: Tetrad
    ) -> dict[Residue3D, Residue3D]:
        """Build the cyclic successor map around a tetrad's base-pairing edges.

        Each base pair connects two edges (Watson-Crick, Hoogsteen, Sugar).
        The ``BasePair3D.score_table`` provides a global preference ordering:
        the direction with the lower score is the "forward" direction
        (donor -> acceptor), matching the established ``cWH``/``cHW`` convention
        and extending it to every Leontis/Westhof pair type.  When a pair is a
        tie (both directions score equally, e.g. ``cWW`` or ``tHH``), its
        direction is left undecided and later inferred by completing the single
        4-cycle from the directed pairs.
        """
        pairs = [tetrad.pair_12, tetrad.pair_23, tetrad.pair_34, tetrad.pair_41]
        successor: dict[Residue3D, Residue3D] = {}
        tie_pairs: list[tuple[Residue3D, Residue3D]] = []

        for pair in pairs:
            score_org = BasePair3D.score_table.get(pair.lw, 20)
            score_rev = BasePair3D.score_table.get(pair.lw.reverse, 20)
            if score_org < score_rev:
                successor[pair.nt1_3d] = pair.nt2_3d
            elif score_rev < score_org:
                successor[pair.nt2_3d] = pair.nt1_3d
            else:
                tie_pairs.append((pair.nt1_3d, pair.nt2_3d))

        if not tie_pairs or len(successor) + len(tie_pairs) != 4:
            return successor

        all_nts = set(tetrad.nucleotides)

        for choices in itertools.product((0, 1), repeat=len(tie_pairs)):
            candidate = dict(successor)
            valid = True
            for (nt1, nt2), choice in zip(tie_pairs, choices, strict=True):
                src, dst = (nt1, nt2) if choice == 0 else (nt2, nt1)
                if src in candidate:
                    valid = False
                    break
                candidate[src] = dst
            if not valid or len(candidate) != 4:
                continue
            if set(candidate.keys()) != all_nts or set(candidate.values()) != all_nts:
                continue
            start = tetrad.nucleotides[0]
            current = start
            for _ in range(4):
                current = candidate.get(current)
                if current is None:
                    break
            if current == start:
                return candidate

        return successor

    def __detect_loop_type(
        self, nt_first: Residue3D, nt_last: Residue3D
    ) -> LoopType | None:
        """Classify a loop by the column transition between its flanking tetrads.

        Columns are ElTetrado tract indices, numbered anticlockwise starting
        from the tetrad-A anchor. A loop is ``+`` when the 3'-flanking column
        is the clockwise-next of the 5'-flanking column (i.e. the column number
        decreases by 1 modulo 4), ``-`` when it is the anticlockwise-next
        (column number increases by 1 modulo 4), and ``diagonal`` (no sign)
        when the two columns are opposite each other.
        """
        column_first = self.__find_tract_index_with_nt(nt_first)
        column_last = self.__find_tract_index_with_nt(nt_last)
        if column_first is None or column_last is None:
            logger.warning(
                "Failed to classify the loop between %s and %s",
                nt_first,
                nt_last,
            )
            return None

        delta = (column_last - column_first) % 4
        if delta == 2:
            return LoopType.diagonal

        sign = "+" if delta == 3 else "-"

        tetrad_with_first = self.__find_tetrad_with_nt(nt_first)
        tetrad_with_last = self.__find_tetrad_with_nt(nt_last)
        if tetrad_with_first is None or tetrad_with_last is None:
            logger.warning(
                "Failed to classify the loop between %s and %s",
                nt_first,
                nt_last,
            )
            return None

        if tetrad_with_first == tetrad_with_last:
            return LoopType.from_value(f"lateral{sign}")
        return LoopType.from_value(f"propeller{sign}")

    # ------------------------------------------------------------------
    # Bulges
    # ------------------------------------------------------------------
    def __find_bulges(self) -> list[BulgeGroup]:
        """
        Detect nucleotides bulging out between two stacked nucleotides of the
        same tract.  Example:
            …-DG13-DC14-DG15-…
                        ↑
                      bulge
        The two flanking nucleotides are chain-consecutive and occupy
        adjacent positions within their tract (genuinely stacked); whatever
        sits between them is the bulge, of any length.
        """
        tract_positions = [
            {nt: position for position, nt in enumerate(tract.nucleotides)}
            for tract in self.tracts
        ]
        bulges: list[BulgeGroup] = []

        for nt_prev, nt_cur, in_between in self.__chain_connectors():
            tract_prev = self.__find_tract_index_with_nt(nt_prev)
            if tract_prev is None or tract_prev != self.__find_tract_index_with_nt(
                nt_cur
            ):
                continue
            if (
                abs(
                    tract_positions[tract_prev][nt_prev]
                    - tract_positions[tract_prev][nt_cur]
                )
                == 1
                and in_between
            ):
                bulges.append(
                    BulgeGroup(
                        in_between,
                        self.__missing_residue_names_between(nt_prev, nt_cur),
                    )
                )

        # keep order, remove duplicate groups with identical nucleotides
        seen: set[frozenset[Residue3D]] = set()
        unique_bulges: list[BulgeGroup] = []
        for bulge in bulges:
            key = frozenset(bulge.nucleotides)
            if key in seen:
                continue
            seen.add(key)
            unique_bulges.append(bulge)
        return unique_bulges

    def __find_snapbacks(self) -> list[Snapback]:
        """
        Detect snapbacks: chain re-entries into an already visited tract.

        Two kinds are recognized:

        - direct reversals: chain-consecutive tetrad nucleotides of the same
          tract whose positions within the tract are not adjacent (e.g. the
          0-nt snapbacks of 5de5/2la7 or the 2-nt one of 2m53); the
          orientation says which side of the connector is the isolated one:
          3' when the re-entering (post-connector) stacked run is shorter -
          its first nucleotide is "injected" into the column - and 5' for the
          mirror case; ties default to 3';
        - re-entries: the chain leaves a tract and returns to it only after
          visiting other columns (the run break is an excursion, not a
          bulge); one of the two adjacent runs is a single nucleotide - the
          "snapped" guanine (e.g. 2mgn G8,G9 -> G24 = 3', 6jcd G6 -> G19 =
          3', 2kpr G1,G2 -> G11 = 5', 6xt7 G1 <- G10,G11,G12 = 5'). The
          orientation follows the chain-terminal run: 5' when the run holding
          the chain's first tetrad nucleotide is the snapped one, 3' when it
          is the run holding the chain's last tetrad nucleotide; ties default
          to the shorter-run rule above. Two stacked stems (two-block
          re-entries such as 2ms9 or the 6jje family) are not snapbacks.
        """
        if len(self.tetrads) == 1:
            return []

        tetrad_nucleotides = sorted(
            [nt for tetrad in self.tetrads for nt in tetrad.nucleotides],
            key=lambda nt: self.global_index[nt],
        )
        tract_positions = [
            {nt: position for position, nt in enumerate(tract.nucleotides)}
            for tract in self.tracts
        ]
        snapbacks: list[Snapback] = []

        for i in range(1, len(tetrad_nucleotides)):
            nt_prev = tetrad_nucleotides[i - 1]
            nt_cur = tetrad_nucleotides[i]
            if nt_prev.chain != nt_cur.chain:
                continue
            tract_index = self.__find_tract_index_with_nt(nt_prev)
            if tract_index is None or tract_index != self.__find_tract_index_with_nt(
                nt_cur
            ):
                continue
            positions = tract_positions[tract_index]
            if abs(positions[nt_prev] - positions[nt_cur]) == 1:
                continue  # stacked — a stack or a bulge connector

            def run_length(
                start: int,
                step: int,
                tract_index: int = tract_index,
                positions: dict = positions,
            ) -> int:
                length = 1
                j = start
                while 0 <= j + step < len(tetrad_nucleotides):
                    nt_a = tetrad_nucleotides[j]
                    nt_b = tetrad_nucleotides[j + step]
                    if nt_a.chain != nt_b.chain:
                        break
                    if self.__find_tract_index_with_nt(nt_b) != tract_index:
                        break
                    if abs(positions[nt_a] - positions[nt_b]) != 1:
                        break
                    length += 1
                    j += step
                return length

            run_prev = run_length(i - 1, -1)
            run_next = run_length(i, 1)
            orientation = (
                SnapbackOrientation.FIVE_PRIME
                if run_next > run_prev
                else SnapbackOrientation.THREE_PRIME
            )
            in_between = list(
                filter(
                    lambda nt: (
                        nt.is_nucleotide
                        and self.global_index[nt_prev]
                        < self.global_index[nt]
                        < self.global_index[nt_cur]
                    ),
                    self.structure3d.residues,
                )
            )
            snapbacks.append(
                Snapback(
                    in_between,
                    orientation,
                    self.__missing_residue_names_between(nt_prev, nt_cur),
                )
            )

        snapbacks.extend(self.__find_re_entry_snapbacks(tetrad_nucleotides))
        return snapbacks

    def __find_re_entry_snapbacks(
        self, tetrad_nucleotides: list[Residue3D]
    ) -> list[Snapback]:
        """Snapbacks formed by tract re-entries after an excursion through
        other columns (see ``__find_snapbacks``)."""
        first_of_chain: dict[str, Residue3D] = {}
        last_of_chain: dict[str, Residue3D] = {}
        for nt in tetrad_nucleotides:
            first_of_chain.setdefault(nt.chain, nt)
            last_of_chain[nt.chain] = nt

        def in_between_nucleotides(nt_prev: Residue3D, nt_cur: Residue3D):
            return list(
                filter(
                    lambda nt: (
                        nt.is_nucleotide
                        and self.global_index[nt_prev]
                        < self.global_index[nt]
                        < self.global_index[nt_cur]
                    ),
                    self.structure3d.residues,
                )
            )

        snapbacks: list[Snapback] = []
        for tract_runs in quadruplex_tract_runs(self):
            for run_prev, run_next in itertools.pairwise(tract_runs):
                chain = run_prev[-1].chain
                if chain != run_next[0].chain:
                    continue  # a chain switch, not a re-entry
                if len(run_prev) > 1 and len(run_next) > 1:
                    continue  # two stacked stems — a two-block, not a snapback

                orientation = None
                for nt in run_prev + run_next:
                    if nt is first_of_chain.get(chain):
                        orientation = SnapbackOrientation.FIVE_PRIME
                        break
                    if nt is last_of_chain.get(chain):
                        orientation = SnapbackOrientation.THREE_PRIME
                        break
                if orientation is None:
                    orientation = (
                        SnapbackOrientation.FIVE_PRIME
                        if len(run_next) > len(run_prev)
                        else SnapbackOrientation.THREE_PRIME
                    )
                snapbacks.append(
                    Snapback(
                        in_between_nucleotides(run_prev[-1], run_next[0]),
                        orientation,
                        self.__missing_residue_names_between(run_prev[-1], run_next[0]),
                    )
                )
        return snapbacks

    def __find_tags(self) -> list[str]:
        """Summarize the quadruplex's structural features as short tags,
        exposed in the text report and the JSON output."""
        if not any(is_pure_g_tetrad(tetrad) for tetrad in self.tetrads):
            return ["non_g4"]

        tags: list[str] = []
        if self.bulges:
            tags.append("has_bulge")
        if self.snapbacks:
            tags.append("has_snapback")
        if self.__has_inverted_strand_polarity():
            tags.append("has_inverted_strand_polarity")
        if self.__has_v_loop():
            tags.append("v_loop")
        if self.__is_two_block():
            tags.append("two_block")
        if self.__is_interlocked():
            tags.append("interlocked")
        if self.loop_class is not None:
            topology_name = TOPOLOGY_NAMES.get(self.loop_class.value)
            if topology_name is not None:
                tags.append(topology_name)
        if self.__has_quadruplex_duplex():
            tags.append("quadruplex_duplex")
        return tags

    def __tetrad_chains(self) -> set[str]:
        return {nt.chain for tetrad in self.tetrads for nt in tetrad.nucleotides}

    def __has_v_loop(self) -> bool:
        """True when the quadruplex contains a V-loop: a cross-column
        connector whose flanking nucleotides lie in different tetrads (a
        propeller-type position) and where one flanking nucleotide is an
        isolated single-nucleotide run of a broken column - the chain snaps
        out of / into that lone guanine (e.g. 2kpr G11-G12, 5zev G20-G21,
        1jjp G10-G11, 2mb4 G3-G5). The loop itself may hold nucleotides
        (the 1-nt V-loops of 2mb4 or the 3qxr family).

        Excluded: same-tetrad 0-nt connectors (lateral/diagonal positions,
        e.g. 4p1d), connectors spanning unmodelled residues (missing loop
        nucleotides, e.g. 7sxp), mirror-image L-nucleotide structures
        (4wb2/4wb3) and non-G4 quadruplexes.
        """
        if self.__is_mirror_image():
            return False
        isolated = isolated_run_nucleotides(quadruplex_tract_runs(self))
        if not isolated:
            return False

        tetrad_of = {
            nt: index
            for index, tetrad in enumerate(self.tetrads)
            for nt in tetrad.nucleotides
        }
        for nt_prev, nt_cur, _ in self.__chain_connectors():
            tract_prev = self.__find_tract_index_with_nt(nt_prev)
            tract_cur = self.__find_tract_index_with_nt(nt_cur)
            if tract_prev is None or tract_cur is None or tract_prev == tract_cur:
                continue
            if tetrad_of[nt_prev] == tetrad_of[nt_cur]:
                continue
            if nt_prev not in isolated and nt_cur not in isolated:
                continue
            if not self.__chain_numbering_is_continuous(nt_prev, nt_cur):
                continue
            return True
        return False

    def __chain_numbering_is_continuous(
        self, nt_prev: Residue3D, nt_cur: Residue3D
    ) -> bool:
        """Residue numbering consistent with the modelled content - excludes
        connectors that span unmodelled (missing) residues."""
        return self.__missing_residue_count_between(nt_prev, nt_cur) == 0

    def __missing_residue_count_between(
        self, nt_prev: Residue3D, nt_cur: Residue3D
    ) -> int:
        """Number of unmodelled residues between two chain-consecutive
        nucleotides: numbering slots minus the residues actually present."""
        index_prev = self.global_index[nt_prev]
        index_cur = self.global_index[nt_cur]
        modelled = sum(
            1
            for residue in self.structure3d.residues
            if residue.chain == nt_prev.chain
            and index_prev < self.global_index[residue] < index_cur
        )
        return abs(nt_cur.number - nt_prev.number) - 1 - modelled

    def __missing_residue_names_between(
        self, nt_prev: Residue3D, nt_cur: Residue3D
    ) -> list[str]:
        """Full names of the residues missing between two chain-consecutive
        nucleotides.

        When the file header reports the unobserved residues (rnapolis >=
        0.16.15) the real identifiers are returned; otherwise placeholders
        (``A.?7`` style) are synthesized from the numbering so that the list
        length always equals the number of missing residues.
        """
        low = min(nt_prev.number, nt_cur.number)
        high = max(nt_prev.number, nt_cur.number)
        slots = range(low + 1, high)
        header_missing = getattr(self.structure3d, "missing_residues", None) or []
        names = [
            Residue(missing.label, missing.auth).full_name
            or f"{missing.chain}.{missing.name}{missing.number}"
            for missing in sorted(header_missing)
            if missing.chain == nt_prev.chain
            and missing.number is not None
            and missing.number in slots
        ]
        if names or not slots:
            return names
        modelled_numbers = {
            residue.number
            for residue in self.structure3d.residues
            if residue.chain == nt_prev.chain and low < residue.number < high
        }
        return [
            f"{nt_prev.chain}.?{number}"
            for number in slots
            if number not in modelled_numbers
        ]

    def __is_mirror_image(self) -> bool:
        """True when all tetrad nucleotides are L-configuration sugars (PDB
        chemical component IDs prefixed with 0, e.g. 0G/0DA in mirror-image
        L-RNA/L-DNA structures such as 4wb2/4wb3)."""
        return all(
            nt.name.startswith("0")
            for tetrad in self.tetrads
            for nt in tetrad.nucleotides
        )

    def __is_interlocked(self) -> bool:
        """True for a multi-chain quadruplex whose chains share ALL columns
        (tracts) and whose arrangement has no clean row partition - the
        columns interpenetrate in a staggered fashion, e.g. the
        interlocked/intercalated dimers such as 1jjp or the locked
        tetramolecular quadruplex pairs such as 1nyd and 1xce. Mutually
        exclusive with two_block, which requires a row partition."""
        if len(self.__tetrad_chains()) < 2:
            return False
        mixed_tracts = sum(
            1
            for tract in self.tracts
            if len({nt.chain for nt in tract.nucleotides}) > 1
        )
        return mixed_tracts == len(self.tracts) and not self.__has_row_partition()

    def __has_inverted_strand_polarity(self) -> bool:
        for tract_polarities in self.strand_polarities:
            if not tract_polarities:
                continue
            has_plus = any(p == StrandPolarity.PLUS for p in tract_polarities)
            has_minus = any(p == StrandPolarity.MINUS for p in tract_polarities)
            if has_plus and has_minus:
                return True
        return False

    def __is_two_block(self) -> bool:
        """True for a quadruplex built from two separate blocks: some column
        is visited by two or more chain segments that each stack at least two
        guanines (the stems share columns), and the blocks stack cleanly - a
        row boundary exists at which every column splits into
        chain-homogeneous halves (unimolecular quadruplexes satisfy this
        vacuously; staggered interlocked arrangements do not). Applies to any
        molecularity (e.g. the bimolecular 6jje/6jjf family and 6sx6/6syk);
        1+2 or 2+1 tetrad caps (Mango, beetroot, sc1 aptamers) are NOT
        two-blocks. 1oz8 is an accepted miss: its real 2+2 arrangement is
        seen as 3 tetrads due to structure quality."""
        if len(self.tetrads) < 2:
            return False

        return self.__has_shared_column_blocks() and self.__has_row_partition()

    def __has_row_partition(self) -> bool:
        """A boundary between stack-adjacent tetrads at which every column
        splits into chain-homogeneous halves exists."""
        return bool(quadruplex_row_partition_boundaries(self))

    def __has_shared_column_blocks(self) -> bool:
        """Vertical split: two chain segments of at least two tetrad
        nucleotides visit the same column."""
        tetrad_nucleotides = sorted(
            [nt for tetrad in self.tetrads for nt in tetrad.nucleotides],
            key=lambda nt: self.global_index[nt],
        )
        tract_positions = [
            {nt: position for position, nt in enumerate(tract.nucleotides)}
            for tract in self.tracts
        ]
        segments: list[list[int]] = [[] for _ in self.tracts]
        run_tract: int | None = None
        run_length = 0
        run_last: Residue3D | None = None

        for nt in tetrad_nucleotides:
            tract_index = self.__find_tract_index_with_nt(nt)
            continues = (
                tract_index is not None
                and tract_index == run_tract
                and run_last is not None
                and run_last.chain == nt.chain
                and abs(
                    tract_positions[tract_index][nt]
                    - tract_positions[tract_index][run_last]
                )
                == 1
            )
            if continues:
                run_length += 1
            else:
                if run_tract is not None:
                    segments[run_tract].append(run_length)
                run_tract = tract_index
                run_length = 1 if tract_index is not None else 0
            run_last = nt
        if run_tract is not None:
            segments[run_tract].append(run_length)

        for tract_segments in segments:
            multi_g_segments = [length for length in tract_segments if length >= 2]
            if len(multi_g_segments) >= 2:
                return True
        return False

    def __has_quadruplex_duplex(self) -> bool:
        """True when at least two canonical Watson-Crick base pairs are formed
        by the quadruplex's loop or bulge nucleotides - a duplex stem growing
        out of a loop (e.g. 6fc9)."""
        loop_or_bulge = {nt for loop in self.loops for nt in loop.nucleotides} | {
            nt for bulge in self.bulges for nt in bulge.nucleotides
        }
        if not loop_or_bulge:
            return False

        seen_pairs = set()
        canonical_pairs = 0
        for base_pair in self.base_pairs:
            if base_pair.lw != LeontisWesthof.cWW:
                continue
            nt1, nt2 = base_pair.nt1_3d, base_pair.nt2_3d
            if nt1 not in loop_or_bulge or nt2 not in loop_or_bulge:
                continue
            pair = frozenset((nt1, nt2))
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)
            pair_letters = {
                nt1.one_letter_name.upper(),
                nt2.one_letter_name.upper(),
            }
            if pair_letters in ({"G", "C"}, {"A", "U"}, {"A", "T"}):
                canonical_pairs += 1
                if canonical_pairs >= 2:
                    return True
        return False

    def __find_tetrad_with_nt(self, nt: Residue3D) -> Tetrad | None:
        for tetrad in self.tetrads:
            if nt in tetrad.nucleotides:
                return tetrad
        return None

    def __find_tract_with_nt(self, nt: Residue3D) -> Tract | None:
        for tract in self.tracts:
            if nt in tract.nucleotides:
                return tract
        return None

    def __find_tract_index_with_nt(self, nt: Residue3D) -> int | None:
        for index, tract in enumerate(self.tracts):
            if nt in tract.nucleotides:
                return index
        return None

    def __classify_by_loops(self) -> LoopClassification | None:
        if len(self.loops) != 3 or any(loop.loop_type is None for loop in self.loops):
            return None

        loop_classes = {
            "ppp": "1",
            "ppl": "2",
            "plp": "3",
            "lpp": "4",
            "pdp": "5",
            "lll": "6",
            "llp": "7",
            "lpl": "8",
            "pll": "9",
            "pdl": "10",
            "ldl": "11",
            "dpd": "12",
            "ldp": "13",
        }
        fingerprint = "".join([loop.loop_type.value[0] for loop in self.loops])  # type: ignore
        if fingerprint not in loop_classes:
            logger.error("Unknown loop classification: %s", fingerprint)
            return None
        subtype = (
            "a"
            if self.loops[0 if fingerprint != "dpd" else 1].loop_type.value[-1] == "-"  # type: ignore
            else "b"
        )
        return LoopClassification.from_value(f"{loop_classes[fingerprint]}{subtype}")

    def __str__(self):
        builder = ""
        if len(self.tetrads) == 1:
            builder += "  single tetrad\n"
            builder += self.__topology_details_str()
            builder += str(self.tetrads[0])
            builder += "    Path:\n"
            builder += f"      {', '.join(self.path)}\n"
        else:
            builder += f"  {self.onzm.value if self.onzm is not None else 'R'}"
            builder += f" {','.join(gba.value for gba in self.gba_classes)}"
            if self.loop_class:
                builder += (
                    f" {self.loop_class.value} {self.loop_class.loop_progression()}"
                )
            else:
                builder += " n/a"
            builder += f" quadruplex with {len(self.tetrads)} tetrads\n"
            if self.tags:
                builder += f"    Tags: {', '.join(self.tags)}\n"
            builder += self.__topology_details_str()
            builder += str(self.tetrad_pairs[0].tetrad1)
            for tetrad_pair in self.tetrad_pairs:
                builder += str(tetrad_pair)
                builder += str(tetrad_pair.tetrad2)
            if self.tracts:
                builder += "\n    Tracts:\n"
                for tract in self.tracts:
                    builder += f"{tract}\n"
            if self.strand_polarities:
                builder += "\n    Strand polarities:\n"
                for tract_index, tract_polarities in enumerate(self.strand_polarities):
                    if not tract_polarities:
                        continue
                    parts = []
                    has_plus = any(p == StrandPolarity.PLUS for p in tract_polarities)
                    has_minus = any(p == StrandPolarity.MINUS for p in tract_polarities)
                    for nt, polarity in zip(
                        self.tracts[tract_index].nucleotides,
                        tract_polarities,
                        strict=True,
                    ):
                        suffix = ""
                        if polarity == StrandPolarity.PLUS:
                            suffix = "+"
                        elif polarity == StrandPolarity.MINUS:
                            suffix = "-"
                        parts.append(f"{nt.full_name}{suffix}")
                    line = "      " + ", ".join(parts)
                    if has_plus and has_minus:
                        line += " (inverted strand polarity)"
                    builder += line + "\n"
            if self.path:
                builder += "\n    Path:\n"
                builder += f"      {', '.join(self.path)}\n"
            if self.loops:
                builder += "\n    Loops:\n"
                for loop in self.loops:
                    builder += f"{loop}\n"
            if self.snapbacks:
                builder += "\n    Snapbacks:\n"
                for snapback in self.snapbacks:
                    builder += f"{snapback}\n"
            if self.bulges:
                builder += "\n    Bulges:\n"
                for bulge in self.bulges:
                    builder += f"      {bulge}\n"
            builder += "\n"
        return builder

    def __topology_details_str(self) -> str:
        builder = ""

        if self.handedness is not None:
            builder += f"    Handedness: {self.handedness.value}\n"
        if self.tetrad_polarities:
            builder += "    Tetrad polarities: "
            builder += ", ".join(
                polarity.value if polarity is not None else "n/a"
                for polarity in self.tetrad_polarities
            )
            builder += "\n"

        return builder


@dataclass
class Helix:
    tetrads: list[Tetrad]
    tetrad_pairs: list[TetradPair]
    structure3d: Structure3D
    global_index: dict[Residue3D, int]
    base_pairs: list[BasePair3D] = field(default_factory=list)
    quadruplexes: list[Quadruplex] = field(init=False)

    def __post_init__(self):
        self.quadruplexes = self.__find_quadruplexes()

    def __find_quadruplexes(self):
        if len(self.tetrad_pairs) == 0:
            return [
                Quadruplex(
                    self.tetrads,
                    [],
                    self.structure3d,
                    self.global_index,
                    self.base_pairs,
                )
            ]

        quadruplexes = []
        tetrads = []

        for tetrad in [self.tetrad_pairs[0].tetrad1] + [
            tetrad_pair.tetrad2 for tetrad_pair in self.tetrad_pairs
        ]:
            if tetrads and tetrad.chains().isdisjoint(tetrads[-1].chains()):
                quadruplexes.append(
                    Quadruplex(
                        tetrads,
                        self.__filter_tetrad_pairs(tetrads),
                        self.structure3d,
                        self.global_index,
                        self.base_pairs,
                    )
                )
                tetrads = []
            tetrads.append(tetrad)

        quadruplexes.append(
            Quadruplex(
                tetrads,
                self.__filter_tetrad_pairs(tetrads),
                self.structure3d,
                self.global_index,
                self.base_pairs,
            )
        )

        return quadruplexes

    def __filter_tetrad_pairs(self, tetrads: list[Tetrad]) -> list[TetradPair]:
        tetrad_set = set(tetrads)
        return [
            tp
            for tp in self.tetrad_pairs
            if tp.tetrad1 in tetrad_set and tp.tetrad2 in tetrad_set
        ]

    def __str__(self):
        builder = ""
        if len(self.tetrads) > 1:
            builder += f"n4-helix with {len(self.tetrads)} tetrads\n"
            for quadruplex in self.quadruplexes:
                builder += str(quadruplex)
        elif len(self.tetrads) == 1:
            builder += "single tetrad without stacking\n"
            builder += str(self.tetrads[0])
        return builder


@dataclass
class Analysis:
    base_interactions: BaseInteractions
    structure3d: Structure3D
    no_reorder: bool
    keep_spurious: bool = False
    global_index: dict[Residue3D, int] = field(init=False)
    mapping: Mapping2D3D = field(init=False)
    tetrads: list[Tetrad] = field(init=False)
    tetrad_scores: dict[Tetrad, dict[Tetrad, "TetradScore"]] = field(init=False)
    tetrad_pairs: list[TetradPair] = field(init=False)
    helices: list[Helix] = field(init=False)
    ions: list[Atom] = field(init=False)
    sugar_puckers: dict[Residue3D, SugarPucker | None] = field(init=False)
    sequence_entries: list[tuple[str, Residue3D | None]] = field(init=False)
    sequence: str = field(init=False)
    line1: str = field(init=False)
    line2: str = field(init=False)

    def __post_init__(self):
        self.global_index = self.__prepare_global_index()
        self.mapping = Mapping2D3D(
            self.structure3d,
            self.base_interactions.base_pairs,
            self.base_interactions.stackings,
            True,
        )
        self.sequence_entries = self.__prepare_sequence_entries()
        self.sequence_entries = self.__prepare_sequence_entries()
        self.tetrads = self.__find_tetrads(self.no_reorder)
        self.tetrad_scores = self.__calculate_tetrad_scores()
        self.tetrad_pairs = self.__find_tetrad_pairs()
        self.helices = self.__find_helices()
        self.sugar_puckers = {
            residue: classify_sugar_pucker(residue)
            for residue in self.structure3d.residues
            if residue.is_nucleotide
        }

        if not self.no_reorder:
            self.__find_best_chain_order()

        if not self.keep_spurious:
            self.helices = self.__filter_spurious_quadruplexes()

        (
            self.sequence,
            self.line1,
            self.line2,
        ) = self.__generate_twoline_dotbracket()
        self.ions = self.__find_ions()
        self.__assign_ions_to_tetrads()

    def __filter_spurious_quadruplexes(self) -> list[Helix]:
        """Discard quadruplexes that fail the spurious quadruplex filter.

        Helices left without a single quadruplex are dropped as well; this
        also removes helices made of a single unstacked tetrad.
        """
        helices = []
        for helix in self.helices:
            quadruplexes = []
            for quadruplex in helix.quadruplexes:
                reasons = quadruplex_spurious_reasons(quadruplex)
                if reasons:
                    logger.info(
                        "Discarding spurious quadruplex (%s tetrads): %s",
                        len(quadruplex.tetrads),
                        "; ".join(reasons),
                    )
                else:
                    quadruplexes.append(quadruplex)
            if quadruplexes:
                helix.quadruplexes = quadruplexes
                helices.append(helix)
        return helices

    def __prepare_sequence_entries(self) -> list[tuple[str, Residue3D | None]]:
        """Full per-chain sequence entries (chain, letter, residue or None)
        for the report and exports, including residues missing from the 3D
        coordinates."""
        return sequence_entries_with_missing(
            sorted(
                filter(lambda nt: nt.is_nucleotide, self.structure3d.residues),
                key=lambda nt: self.global_index[nt],
            ),
            getattr(self.structure3d, "missing_residues", None),
            getattr(self.structure3d, "sequence_by_entity", None),
        )

    def __prepare_global_index(self) -> dict[Residue3D, int]:
        result = {}
        i = 1
        for residue in self.structure3d.residues:
            result[residue] = i
            i += 1
        return result

    def __find_tetrads(self, no_reorder=False) -> list[Tetrad]:
        # search for a tetrad: i -> j -> k -> l
        #                      ^--------------^
        tetrads = []
        for i in self.mapping.base_pair_graph:
            for j in self.mapping.base_pair_graph[i]:
                if j == i:
                    continue
                for k in self.mapping.base_pair_graph[j]:
                    if k in (i, j):
                        continue
                    for n4 in self.mapping.base_pair_graph[k]:
                        if n4 in (i, j, k) or i not in self.mapping.base_pair_graph[n4]:
                            continue
                        if Tetrad.is_valid(i, j, k, n4, self.mapping.base_pair_dict):
                            pair_12 = self.mapping.base_pair_dict[(i, j)]
                            pair_23 = self.mapping.base_pair_dict[(j, k)]
                            pair_34 = self.mapping.base_pair_dict[(k, n4)]
                            pair_41 = self.mapping.base_pair_dict[(n4, i)]
                            tetrads.append(
                                Tetrad(
                                    i,
                                    j,
                                    k,
                                    n4,
                                    pair_12,
                                    pair_23,
                                    pair_34,
                                    pair_41,
                                    self.global_index,
                                )
                            )

        # build graph of tetrads
        while tetrads:
            graph = defaultdict(list)
            for ti, tj in itertools.combinations(tetrads, 2):
                if not ti.is_disjoint(tj):
                    graph[ti].append(tj)
                    graph[tj].append(ti)

            # remove tetrad which conflicts the most with others
            # in case of a tie, remove one which has the worst planarity deviation (highest rmsd)
            candidates = sorted(
                tetrads,
                key=lambda t: (len(graph[t]), t.planarity_deviation["rmsd"]),
                reverse=True,
            )
            if len(graph[candidates[0]]) > 0:
                tetrads.remove(candidates[0])
            else:
                break

        return sorted(
            tetrads,
            key=lambda t: min(self.global_index[nt] for nt in t.nucleotides),
        )

    def __calculate_tetrad_scores(
        self,
    ) -> dict[Tetrad, dict[Tetrad, tuple[int, tuple, tuple]]]:
        """Score every tetrad pair to decide how tetrads should be connected.

        Each candidate correspondence combines two independent signals:
        sequential proximity in 5' order and observed stacking interactions. The
        winning permutation provides the residue-to-residue mapping later used to
        build tracts and adjacent ``TetradPair`` objects.
        """

        def is_next_by_stacking(nt1: Residue3D, nt2: Residue3D) -> bool:
            return nt2 in self.mapping.stacking_graph.get(nt1, [])

        def is_next_sequentially(nt1: Residue3D, nt2: Residue3D) -> float:
            """
            Return a *fractional* sequential-proximity score between two residues.

            1.0  - same chain, immediately adjacent (|idx₁ - idx₂| == 1)
            0.5  - same chain, one residue apart  (|idx₁ - idx₂| == 2)
            0.0  - anything else
            """
            if nt1.chain != nt2.chain:
                return 0.0

            # both residues must be indexed; otherwise return 0
            idx1 = self.global_index.get(nt1)
            idx2 = self.global_index.get(nt2)
            if idx1 is None or idx2 is None:
                logger.warning(
                    "Residue missing in global_index while scoring sequential proximity"
                )
                return 0.0

            diff = abs(idx1 - idx2)
            if diff == 1:
                return 1.0
            if diff == 2:
                return 0.5
            return 0.0

        tetrad_scores: dict[
            Tetrad,
            dict[Tetrad, tuple[int, tuple[Residue3D, ...], tuple[Residue3D, ...]]],
        ] = defaultdict(dict)

        for ti, tj in itertools.combinations(self.tetrads, 2):
            best_score_sequential = 0
            best_score_stacking = 0

            nts1 = ti.nucleotides
            nts2 = tj.nucleotides
            n1, n2, n3, n4 = nts2
            viable_permutations = [
                (n1, n2, n3, n4),
                (n2, n3, n4, n1),
                (n3, n4, n1, n2),
                (n4, n1, n2, n3),
                (n1, n4, n3, n2),
                (n4, n3, n2, n1),
                (n3, n2, n1, n4),
                (n2, n1, n4, n3),
            ]

            for permutation in viable_permutations:
                flags_sequential: list[float] = [
                    is_next_sequentially(nts1[i], permutation[i]) for i in range(4)
                ]
                flags_stacking: list[int] = [
                    int(is_next_by_stacking(nts1[i], permutation[i])) for i in range(4)
                ]

                # compute independent scores
                score_sequential = sum(flags_sequential)
                score_stacking = sum(flags_stacking)

                if (score_sequential, score_stacking) > (
                    best_score_sequential,
                    best_score_stacking,
                ):
                    best_score_sequential, best_score_stacking = (
                        score_sequential,
                        score_stacking,
                    )
                    nts2 = permutation

            tetrad_scores[ti][tj] = TetradScore(
                best_score_sequential,
                best_score_stacking,
                nts1,
                nts2,
            )
            tetrad_scores[tj][ti] = TetradScore(
                best_score_sequential,
                best_score_stacking,
                nts2,
                nts1,
            )

        # log information about tetrad scores
        logger.debug("Tetrad scores:")
        for ti in self.tetrads:
            logger.debug("%r", ti)
            for tj in self.tetrads:
                if ti != tj:
                    logger.debug("\t%r\t%s", tj, tetrad_scores[ti][tj])

        return tetrad_scores

    def __find_tetrad_pairs(self) -> list[TetradPair]:
        """Construct a best adjacent-stack path through the detected tetrads.

        The resulting order is the internal stack order used for tract extension,
        adjacent pair geometry, and helix grouping. It is intentionally separate
        from tetrad 5' labeling and from g4composer build order.
        """

        def next_tetrad_scoring(
            ti: Tetrad, tj: Tetrad, candidates: Iterable[Tetrad]
        ) -> tuple[int, float, int, int, int]:
            """
            Provide a sorting key that expresses how "good" a follow-up tetrad *tj*
            is when the current end of a tentative helix is *ti*.

            The tuple is evaluated lexicographically by ``max(…, key=next_tetrad_scoring)``.
            Order of importance:
              1. stacking alignment component  (0-4)
              2. negative distance between tetrad centroids (favour closer tetrads)
              3. sequential alignment component (0-4)
              4. negative sum of *tj*'s scores to the remaining *candidates*
                 - favours tetrads that are already **less** compatible with the
                   rest of the pool so we can grow one continuous helix first.
              5. deterministic tiebreaker = negative original index of *tj*
            """
            score_direct = self.tetrad_scores[ti].get(tj)

            # direct scores (0-4 each)
            sequential_direct = score_direct.sequential if score_direct else 0
            stacking_direct = score_direct.stacking if score_direct else 0

            # geometric distance between tetrad centroids
            distance = -float(numpy.linalg.norm(tj.center() - ti.center()))

            # how strongly tj is connected to still-unvisited tetrads
            total_candidates = -sum(
                (
                    self.tetrad_scores[tj][tk].sequential
                    + self.tetrad_scores[tj][tk].stacking
                    if tk in self.tetrad_scores[tj]
                    else 0
                )
                for tk in candidates
            )

            return (
                stacking_direct,
                distance,
                sequential_direct,
                total_candidates,
                -self.tetrads.index(tj),
            )

        tetrads = list(self.tetrads)
        best_score = [0, 0]
        best_order = tetrads

        for ti in tetrads:
            score = [0, 0]
            order = [ti]
            candidates = set(self.tetrads) - {ti}
            current = ti

            while candidates:
                tj = max(
                    candidates,
                    key=lambda tk: next_tetrad_scoring(current, tk, candidates),
                )
                score[0] += self.tetrad_scores[current][tj].sequential
                score[1] += self.tetrad_scores[current][tj].stacking
                order.append(tj)
                candidates.remove(tj)
                current = tj

            if score > best_score:
                best_score = score
                best_order = order

        tetrad_pairs = []

        for i in range(1, len(best_order)):
            ti, tj = best_order[i - 1], best_order[i]
            score_seq = self.tetrad_scores[ti][tj].sequential
            score_stack = self.tetrad_scores[ti][tj].stacking

            if (score_seq > 0) or (score_stack > 0):
                nts1, nts2 = (
                    self.tetrad_scores[ti][tj].nts1,
                    self.tetrad_scores[ti][tj].nts2,
                )
                stacked = {nts1[i]: nts2[i] for i in range(4)}
                stacked.update({v: k for k, v in stacked.items()})
                tetrad_pairs.append(TetradPair(ti, tj, stacked, self.global_index))
                stacked_order = (
                    stacked[ti.nt1],
                    stacked[ti.nt2],
                    stacked[ti.nt3],
                    stacked[ti.nt4],
                )
                tj.reorder_to_match_other_tetrad(stacked_order)

        return tetrad_pairs

    def __find_helices(self):
        helices = []
        helix_tetrads = []
        helix_tetrad_pairs = []

        for tp in self.tetrad_pairs:
            ti, tj = tp.tetrad1, tp.tetrad2
            if not helix_tetrads:
                helix_tetrads.append(ti)
            score_seq = self.tetrad_scores[helix_tetrads[-1]][tj].sequential
            score_stack = self.tetrad_scores[helix_tetrads[-1]][tj].stacking
            if (score_seq > 0) or (score_stack > 0):
                helix_tetrads.append(tj)
                helix_tetrad_pairs.append(tp)
            else:
                helices.append(
                    Helix(
                        helix_tetrads,
                        helix_tetrad_pairs,
                        self.structure3d,
                        self.global_index,
                        self.mapping.base_pairs,
                    )
                )
                helix_tetrads = [ti, tj]
                helix_tetrad_pairs = [tp]

        if helix_tetrads:
            helices.append(
                Helix(
                    helix_tetrads,
                    helix_tetrad_pairs,
                    self.structure3d,
                    self.global_index,
                    self.mapping.base_pairs,
                )
            )

        for tetrad in self.tetrads:
            if not any(tetrad in helix.tetrads for helix in helices):
                helices.append(
                    Helix(
                        [tetrad],
                        [],
                        self.structure3d,
                        self.global_index,
                        self.mapping.base_pairs,
                    )
                )

        return helices

    def __find_best_chain_order(self):
        chain_groups = self.__group_related_chains()
        final_order = []

        for chains in chain_groups:
            best_permutation, best_score = tuple(chains), (1e10, 1e10)

            if len(chains) > 1:
                for permutation in itertools.permutations(chains):
                    self.__reorder_chains(permutation)
                    classifications = [t.onz for h in self.helices for t in h.tetrads]
                    logger.debug(
                        "Checking reorder: %s %s",
                        " ".join(permutation),
                        " ".join(c.value for c in classifications),
                    )

                    onz_score = sum(c.score() for c in classifications)
                    chain_order_score = self.__chain_order_score(permutation)
                    score = (onz_score, chain_order_score)

                    if score < best_score:
                        best_score = score
                        best_permutation = permutation
                    elif score == best_score:
                        # in case of a tie, pick permutation earlier in lexicographical sense
                        best_permutation = min(best_permutation, permutation)

            final_order.extend(best_permutation)

        if len(final_order) > 1:
            self.__reorder_chains(final_order)
            classifications = [t.onz for h in self.helices for t in h.tetrads]
            logger.debug(
                "Selected chain order: %s %s",
                " ".join(final_order),
                " ".join(onz.value for onz in classifications),
            )

            self.tetrads = self.__find_tetrads(True)
            self.tetrad_scores = self.__calculate_tetrad_scores()
            self.tetrad_pairs = self.__find_tetrad_pairs()
            self.helices = self.__find_helices()

    def __group_related_chains(self) -> list[list[str]]:
        candidates_prior: set[frozenset[str]] = set()

        for h in self.helices:
            for t in h.tetrads:
                candidates_prior.add(
                    frozenset([t.nt1.chain, t.nt2.chain, t.nt3.chain, t.nt4.chain])
                )

        candidates: list[set[str]] = [set(c) for c in candidates_prior]
        changed = True

        while changed:
            changed = False

            for i, j in itertools.combinations(range(len(candidates)), 2):
                qi, qj = candidates[i], candidates[j]

                if not qi.isdisjoint(qj):
                    qi.update(qj)
                    del candidates[j]
                    changed = True
                    break

        candidates = sorted(candidates, key=len, reverse=True)
        groups: list[set[str]] = []

        for candidate in candidates:
            if any(group.issuperset(candidate) for group in groups):
                continue
            groups.append(candidate)

        return sorted([sorted(group) for group in groups], key=lambda x: x[0])

    def __reorder_chains(self, chain_order: Iterable[str]):
        i = 1
        for chain in chain_order:
            for nt in self.structure3d.residues:
                if nt.chain == chain:
                    self.global_index[nt] = i
                    i += 1
        for nt in self.structure3d.residues:
            if nt.chain not in chain_order:
                self.global_index[nt] = i
                i += 1

        if len(self.tetrad_pairs) > 0:
            self.tetrad_pairs[0].tetrad1.reorder_to_match_5p_3p()
            for tp in self.tetrad_pairs:
                order = (
                    tp.stacked[tp.tetrad1.nt1],
                    tp.stacked[tp.tetrad1.nt2],
                    tp.stacked[tp.tetrad1.nt3],
                    tp.stacked[tp.tetrad1.nt4],
                )
                tp.tetrad2.reorder_to_match_5p_3p()  # this is required to recalculate ONZ
                tp.tetrad2.reorder_to_match_other_tetrad(order)

    def __chain_order_score(self, chain_order: tuple[str, ...]) -> int:
        chain_pairs = []
        for h in self.helices:
            for t in h.tetrads:
                for p in [t.pair_12, t.pair_23, t.pair_34, t.pair_41]:
                    c1 = p.nt1.chain
                    c2 = p.nt2.chain
                    if c1 != c2 and c1 in chain_order and c2 in chain_order:
                        chain_pairs.append([c1, c2])
        sum_sq = 0
        for c1, c2 in chain_pairs:
            sum_sq += (chain_order.index(c1) - chain_order.index(c2)) ** 2
        return sum_sq

    def __find_ions(self) -> list[Atom]:
        metal_atom_names = {ion.value.upper() for ion in Ion}
        ions = []
        used = set()
        for residue in self.structure3d.residues:
            for atom in residue.atoms:
                if atom.name.upper() in metal_atom_names:
                    coordinates = tuple(atom.coordinates)
                    if coordinates not in used:
                        ions.append(atom)
                        used.add(coordinates)
        return ions

    def __assign_ions_to_tetrads(self):
        if len(self.tetrads) == 0:
            return

        ions_channel = defaultdict(list)
        ions_outside = defaultdict(list)

        for ion in self.ions:
            min_distance = math.inf
            min_tetrad = self.tetrads[0]

            for tetrad in self.tetrads:
                distance = numpy.linalg.norm(ion.coordinates - tetrad.center()).item()
                if distance < min_distance:
                    min_distance = distance
                    min_tetrad = tetrad

            # TODO: verify threshold of 6A between an ion and tetrad channel
            if min_distance < 6.0:
                ions_channel[min_tetrad].append(ion)
                continue

            min_distance = math.inf
            min_tetrad = self.tetrads[0]
            min_nt = min_tetrad.nt1

            for tetrad in self.tetrads:
                for nt in tetrad.nucleotides:
                    for atom in nt.atoms:
                        distance = numpy.linalg.norm(
                            ion.coordinates - atom.coordinates
                        ).item()
                        if distance < min_distance:
                            min_distance = distance
                            min_tetrad = tetrad
                            min_nt = nt

            # TODO: verify threshold of 3A between an ion and an atom
            if min_distance < 3.0:
                ions_outside[(min_tetrad, min_nt)].append(ion)
                continue

            logger.debug(
                "Skipping an ion, because it is too far from any tetrad (distance=%s)",
                min_distance,
            )

        for tetrad, ions in ions_channel.items():
            tetrad.ions_channel = ions
        for pair, ions in ions_outside.items():
            tetrad, residue = pair
            tetrad.ions_outside[residue] = ions

    def __generate_twoline_dotbracket(
        self,
    ) -> tuple[str, str, str, dict[Residue3D, int]]:
        layer1, layer2 = [], []

        for tetrad in self.tetrads:
            tetrad_copy = copy.deepcopy(tetrad)
            tetrad_copy.reorder_to_match_5p_3p()
            score_org = BasePair3D.score_table.get(tetrad_copy.pair_12.lw, 20)
            score_rev = BasePair3D.score_table.get(tetrad_copy.pair_12.lw.reverse, 20)

            if score_org < score_rev:
                layer1.extend([tetrad_copy.pair_12, tetrad_copy.pair_34])
                layer2.extend([tetrad_copy.pair_23, tetrad_copy.pair_41])
            else:
                layer1.extend([tetrad_copy.pair_23, tetrad_copy.pair_41])
                layer2.extend([tetrad_copy.pair_12, tetrad_copy.pair_34])

        sequence, line1 = self.__dot_bracket(layer1)
        _, line2 = self.__dot_bracket(layer2)
        return sequence, line1, line2

    def __dot_bracket(self, pairs: list[BasePair3D]) -> tuple[str, str]:
        bpseq_index = {}
        bpseq_entries = []

        for offset, (_, letter, residue) in enumerate(self.sequence_entries):
            if residue is not None:
                bpseq_index[residue] = offset + 1
            bpseq_entries.append(Entry(offset + 1, letter, 0))

        for pair in pairs:
            ni = bpseq_index[pair.nt1_3d]
            nj = bpseq_index[pair.nt2_3d]
            bpseq_entries[ni - 1].pair = nj
            bpseq_entries[nj - 1].pair = ni

        structure_rendering = BpSeq(bpseq_entries).fcfs.structure

        sequence_chars = []
        structure_chars = []
        previous_chain = None
        for offset, (chain, letter, residue) in enumerate(self.sequence_entries):
            if previous_chain is not None and chain != previous_chain:
                sequence_chars.append("-")
                structure_chars.append("-")
            sequence_chars.append(letter)
            structure_chars.append(
                "?" if residue is None else structure_rendering[offset]
            )
            previous_chain = chain
        return "".join(sequence_chars), "".join(structure_chars)

    def __str__(self):
        builder = f"Chain order: {' '.join(self.__chain_order())}\n"
        for helix in self.helices:
            builder += str(helix)
        builder += f"{self.sequence}\n{self.line1}\n{self.line2}"
        return builder

    def __chain_order(self) -> list[str]:
        only_nucleic_acids = filter(
            lambda nt: nt.is_nucleotide, self.structure3d.residues
        )
        return list(
            {
                nt.chain: 0
                for nt in sorted(
                    only_nucleic_acids,
                    key=lambda nt: self.global_index[nt],
                )
            }.keys()
        )

    def __is_conflicted(self, bp1: BasePair3D, bp2: BasePair3D) -> bool:
        xi, yi = sorted(
            [
                self.global_index.get(bp1.nt1_3d, sys.maxsize),
                self.global_index.get(bp1.nt2_3d, sys.maxsize),
            ]
        )
        xj, yj = sorted(
            [
                self.global_index.get(bp2.nt1_3d, sys.maxsize),
                self.global_index.get(bp2.nt2_3d, sys.maxsize),
            ]
        )
        return xi < xj < yi < yj or xj < xi < yj < yi

    def canonical(self) -> list[BasePair3D]:
        return [
            base_pair for base_pair in self.mapping.base_pairs if base_pair.is_canonical
        ]

    def is_basepair_in_tetrad(self, base_pair: BasePair3D) -> bool:
        for tetrad in self.tetrads:
            if base_pair in (
                tetrad.pair_12,
                tetrad.pair_23,
                tetrad.pair_34,
                tetrad.pair_41,
                tetrad.pair_12.reverse,
                tetrad.pair_23.reverse,
                tetrad.pair_34.reverse,
                tetrad.pair_41.reverse,
            ):
                return True
        return False


@dataclass
class Visualizer:
    analysis: Analysis
    tetrads: list[Tetrad]
    complete2d: bool
    global_index: dict[Residue3D, int]
    onz_dict: dict[BasePair3D, ONZ] = field(init=False)

    def __post_init__(self):
        self.onz_dict = {}

    def visualize(self, prefix: str, suffix: str, output_dir: str | None = None):
        # the temporary fasta is consumed by the R subprocess below and removed
        # on process exit, so it deliberately stays open without a context manager
        fasta = tempfile.NamedTemporaryFile("w+", suffix=".fasta")  # noqa: SIM115
        fasta.write(f">{prefix}-{suffix}\n")
        # keep the fasta in the modelled-residue index space used by the .helix
        # file below (the enriched report sequence would shift the coordinates)
        fasta.write("".join(nt.one_letter_name for nt in self.__modelled_nucleotides()))
        fasta.flush()

        layer1, layer2 = [], []
        for tetrad in self.tetrads:
            tetrad_copy = copy.deepcopy(tetrad)
            tetrad_copy.reorder_to_match_5p_3p()
            score_org = BasePair3D.score_table.get(tetrad_copy.pair_12.lw, 20)
            score_rev = BasePair3D.score_table.get(tetrad_copy.pair_12.lw.reverse, 20)

            if score_org < score_rev:
                layer1.extend([tetrad_copy.pair_12, tetrad_copy.pair_34])
                layer2.extend([tetrad_copy.pair_23, tetrad_copy.pair_41])
            else:
                layer1.extend([tetrad_copy.pair_23, tetrad_copy.pair_41])
                layer2.extend([tetrad_copy.pair_12, tetrad_copy.pair_34])

            self.onz_dict[tetrad_copy.pair_12] = tetrad.onz
            self.onz_dict[tetrad_copy.pair_23] = tetrad.onz
            self.onz_dict[tetrad_copy.pair_34] = tetrad.onz
            self.onz_dict[tetrad_copy.pair_41] = tetrad.onz

        helix1 = self.__to_helix(
            layer1, self.analysis.canonical() if self.complete2d else []
        )
        helix2 = self.__to_helix(layer2)

        currdir = os.path.dirname(os.path.realpath(__file__))
        output_pdf = f"{prefix}-{suffix}.pdf"
        if output_dir:
            output_pdf = os.path.join(output_dir, output_pdf)
        run = subprocess.run(
            [
                os.path.join(currdir, "quadraw.R"),
                fasta.name,
                helix1.name,
                helix2.name,
                output_pdf,
            ],
            capture_output=True,
            check=False,
        )
        if run.returncode == 0:
            print("\nPlot:", output_pdf)
        else:
            logger.error(
                "Failed to prepare visualization, reason:\n  %s",
                run.stderr.decode(),
            )

    def __modelled_nucleotides(self) -> list[Residue3D]:
        return [nt for nt in self.analysis.structure3d.residues if nt.is_nucleotide]

    def __to_helix(
        self, layer: list[BasePair3D], canonical: list[BasePair3D] | None = None
    ) -> IO[str]:
        onz_value = {
            ONZ.O_PLUS: 1,
            ONZ.O_MINUS: 2,
            ONZ.N_PLUS: 3,
            ONZ.N_MINUS: 4,
            ONZ.Z_PLUS: 5,
            ONZ.Z_MINUS: 6,
        }
        nucleotides = self.__modelled_nucleotides()
        chain = None
        shifts = {}
        shift_value = 0

        for nucleotide in nucleotides:
            if chain and chain != nucleotide.chain:
                shift_value += 1
            shifts[nucleotide] = shift_value
            chain = nucleotide.chain

        # the .helix file is returned to the caller (consumed by the R
        # subprocess), so it deliberately stays open without a context manager
        helix = tempfile.NamedTemporaryFile("w+", suffix=".helix")  # noqa: SIM115
        helix.write(f"#{len(self.__modelled_nucleotides()) + 1}\n")
        helix.write("i\tj\tlength\tvalue\n")

        for pair in layer:
            x, y = (
                nucleotides.index(pair.nt1_3d) + 1 + shifts[pair.nt1_3d],
                nucleotides.index(pair.nt2_3d) + 1 + shifts[pair.nt2_3d],
            )
            onz = self.onz_dict[pair]
            helix.write(f"{x}\t{y}\t1\t{onz_value.get(onz, 7)}\n")
        if canonical:
            for pair in canonical:
                x, y = (
                    nucleotides.index(pair.nt1_3d) + 1 + shifts[pair.nt1_3d],
                    nucleotides.index(pair.nt2_3d) + 1 + shifts[pair.nt2_3d],
                )
                helix.write(f"{x}\t{y}\t1\t8\n")

        helix.flush()
        return helix


class AnalysisSimple:
    def __init__(self, base_interactions: BaseInteractions, structure3d: Structure3D):
        self.mapping = Mapping2D3D(
            structure3d,
            base_interactions.base_pairs,
            base_interactions.stackings,
            False,
        )

    def has_tetrads(self) -> bool:
        tetrads = set()
        for i in self.mapping.base_pair_graph:
            for j in self.mapping.base_pair_graph[i]:
                if j == i:
                    continue
                for k in self.mapping.base_pair_graph[j]:
                    if k in (i, j):
                        continue
                    for n4 in self.mapping.base_pair_graph[k]:
                        if n4 in (i, j, k) or i not in self.mapping.base_pair_graph[n4]:
                            continue
                        if Tetrad.is_valid(i, j, k, n4, self.mapping.base_pair_dict):
                            tetrads.add(frozenset([i, j, k, n4]))
                        if len(tetrads) > 1:
                            return True
        return False


def center_of_mass(atoms: list[Atom]) -> numpy.typing.NDArray[numpy.floating]:
    coords = [atom.coordinates for atom in atoms]
    xs = (coord[0] for coord in coords)
    ys = (coord[1] for coord in coords)
    zs = (coord[2] for coord in coords)
    return numpy.array(
        (sum(xs) / len(coords), sum(ys) / len(coords), sum(zs) / len(coords))
    )


def eltetrado(
    base_interactions: BaseInteractions,
    structure3d: Structure3D,
    no_reorder: bool,
    keep_spurious: bool = False,
) -> Analysis:
    return Analysis(base_interactions, structure3d, no_reorder, keep_spurious)


def has_tetrad(base_interactions: BaseInteractions, structure3d: Structure3D) -> bool:
    structure = AnalysisSimple(base_interactions, structure3d)
    return structure.has_tetrads()
