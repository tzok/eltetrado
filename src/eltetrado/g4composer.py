"""Generate plain-text exports for g4composer.

This module has to bridge several ordering conventions used elsewhere in ElTetrado:

- ``analysis.global_index`` defines nucleotide 5'->3' order across chains,
- ``quadruplex.tetrads`` / ``quadruplex.tetrad_pairs`` follow adjacent stack order,
- tetrad letters in ``quadruplex.path`` are assigned by tetrad 5'->3' order,
- g4composer builds a quadruplex in the order implied by first appearance of
  tetrad letters in ``path``,
- g4composer path columns are numbered clockwise, whereas ElTetrado tracts are
  numbered anticlockwise after the first column anchor.

The exporter therefore keeps ElTetrado topology intact, but remaps it at export
time into g4composer's build-order and path-numbering conventions.

Two export modes are available (``G4ComposerMode``):

- ``CLASSIC`` (default) keeps physical tract ids as path column numbers (with
  the clockwise remap) and exports twist as a best-fit rotation of all four
  matched tract C1' sets around one shared quadruplex-wide axis.
- ``CONSISTENT`` keeps the classic column numbering but shifts it cyclically
  within each tetrad so that the first tetrad nucleotide visited by the chain
  becomes number 1 (classic ``A1;...;B4;...;B3;...;B2;...;B1`` becomes
  ``A1;...;B1;...;B4;...;B3;...;B2``), and exports twist as the torsion angle
  C1'(1) - center(tetrad 1) - center(tetrad 2) - C1'(1) between consecutive
  build-order tetrads, signed by the interval's rise so it follows the fixed,
  rise-positive quadruplex axis. Rise and every other exported field are
  mode-independent; the shared axis is oriented so the first build-order rise
  is positive in both modes.
"""

import itertools
import logging
import math
import os
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum

import numpy
from rnapolis.common import BpSeq, Entry, GlycosidicBond, Molecule
from rnapolis.tertiary import (
    Residue3D,
    calculate_torsion_angle_coords,
)

from eltetrado.analysis import (
    Analysis,
    Quadruplex,
    Tetrad,
    calculate_best_fit_rotation_around_axis,
    collect_nucleobase_atoms,
    get_plane_vectors,
    residue_c1prime_coordinates,
    sequence_entries_with_missing,
)
from eltetrado.model import SugarPucker

logger = logging.getLogger(__name__)


class G4ComposerMode(Enum):
    """Path/twist numbering convention used by the g4composer exporter.

    See the module docstring for the full description of both modes.
    """

    CLASSIC = "classic"
    CONSISTENT = "consistent"


@dataclass(frozen=True)
class G4ComposerEntry:
    name: str
    sequence: str
    structure: str
    chi: str
    sugar: str
    orient: str
    rise: str
    twist: str
    path: str

    def serialize(self) -> str:
        fields = [
            ("name", self.name),
            ("sequence", self.sequence),
            ("structure", self.structure),
            ("chi", self.chi),
            ("sugar", self.sugar),
            ("orient", self.orient),
            ("rise", self.rise),
            ("twist", self.twist),
            ("path", self.path),
        ]
        return "\n".join(f"{key:<11} {value}" for key, value in fields) + "\n"


def input_name(path: str) -> str:
    basename = os.path.basename(path)
    root, ext = os.path.splitext(basename)
    if ext == ".gz":
        root, _ = os.path.splitext(root)
    return root


def eligible_quadruplexes(analysis: Analysis) -> list[Quadruplex]:
    """Return g4composer-eligible quadruplexes.

    g4composer currently supports a single unimolecular quadruplex only, so we
    filter out quadruplexes built from multiple chains here.
    """
    return [
        quadruplex
        for helix in analysis.helices
        for quadruplex in helix.quadruplexes
        if len(quadruplex.tetrads) >= 2 and len(export_chains(quadruplex)) == 1
    ]


def generate_g4composer_entry(
    analysis: Analysis,
    quadruplex: Quadruplex,
    name: str,
    mode: G4ComposerMode = G4ComposerMode.CLASSIC,
) -> G4ComposerEntry:
    entries = export_entries(analysis, quadruplex)
    build_order = tetrads_in_build_order(quadruplex)
    tetrad_labels = quadruplex.tetrad_labels_by_5p_order()
    tetrad_to_label = {
        tetrad: tetrad_labels[index] for index, tetrad in enumerate(quadruplex.tetrads)
    }

    return G4ComposerEntry(
        name=name,
        sequence=export_sequence(entries, molecule_by_chain(entries)),
        structure=export_structure(analysis, quadruplex, entries),
        chi=export_chi(entries),
        sugar=export_sugar(analysis, entries),
        orient=export_orient(quadruplex, build_order, tetrad_to_label),
        rise=export_rise(quadruplex),
        twist=export_twist(quadruplex, mode),
        path=export_path(quadruplex, mode),
    )


def export_residues(analysis: Analysis, quadruplex: Quadruplex) -> list[Residue3D]:
    """Export full polymer chains participating in the selected quadruplex.

    g4composer expects the quadruplex together with its 5' and 3' flanking
    single-stranded context, so we export every nucleotide from the involved
    chain(s), not only residues that belong to tetrads. Ligand nucleotides
    (GTP, cGMP, ...) are excluded; some of them are re-added by
    ``export_entries`` when they extend the chain termini.
    """
    chains = export_chains(quadruplex)
    return [
        residue
        for residue in sorted(
            analysis.structure3d.residues,
            key=lambda residue: analysis.global_index[residue],
        )
        if residue.chain in chains and residue in analysis.polymer_residues
    ]


def export_chains(quadruplex: Quadruplex) -> list[str]:
    return sorted(
        {
            nt.chain
            for tetrad in quadruplex.tetrads
            for nt in tetrad.nucleotides
            if nt.is_nucleotide
        }
    )


def export_entries(
    analysis: Analysis, quadruplex: Quadruplex
) -> list[tuple[str, str, Residue3D | None, bool]]:
    """Full chain-span entries for the export: modelled polymer residues
    interleaved with residues missing from the coordinates (their real
    one-letter codes when the header sequence is available, else '?').

    Ligand nucleotides that covalently extend the involved chain (5'/3'
    terminal nucleosides, e.g. 8wa4) or form a canonical base pair with its
    terminal residue (e.g. the 7oa3 GTP ligand) are re-added as extra
    sequence letters at the appropriate sequence end."""
    entries = sequence_entries_with_missing(
        export_residues(analysis, quadruplex),
        getattr(analysis.structure3d, "missing_residues", None),
        getattr(analysis.structure3d, "sequence_by_entity", None),
    )
    return extend_entries_with_terminal_ligands(analysis, quadruplex, entries)


def terminal_ligand_letters(
    analysis: Analysis, quadruplex: Quadruplex, entries
) -> list[tuple[int, tuple[str, str, Residue3D | None, bool]]]:
    """Compute ligand entries to insert into the export span.

    Returns (insertion index, entry) pairs. Insertion positions are indices
    into the ligand-free ``entries`` list, so adding them from the largest
    index first keeps earlier positions valid.
    """
    chains = export_chains(quadruplex)
    span_residues = [residue for _, _, residue, _ in entries if residue is not None]
    if not span_residues:
        return []
    first_by_chain: dict[str, Residue3D] = {}
    last_by_chain: dict[str, Residue3D] = {}
    for residue in span_residues:
        first_by_chain.setdefault(residue.chain, residue)
        last_by_chain[residue.chain] = residue

    ligands = [
        residue
        for residue in analysis.structure3d.residues
        if residue.is_nucleotide
        and residue.chain in chains
        and residue not in analysis.polymer_residues
    ]
    if not ligands:
        return []

    canonical = analysis.canonical()
    block_start: dict[str, int] = {}
    block_end: dict[str, int] = {}
    for index, (chain, _, residue, _) in enumerate(entries):
        if residue is not None:
            block_start.setdefault(chain, index)
            block_end[chain] = index + 1

    insertions = []
    for ligand in ligands:
        # covalent extension: the ligand nucleoside bonds to a polymer residue
        bonded = None
        for partner in (
            first_by_chain.get(ligand.chain),
            last_by_chain.get(ligand.chain),
        ):
            if partner is not None and (
                partner.is_connected(ligand) or ligand.is_connected(partner)
            ):
                bonded = partner
                break
        if bonded is None:
            bonded = next(
                (
                    partner
                    for partner in span_residues
                    if partner.chain == ligand.chain
                    and (partner.is_connected(ligand) or ligand.is_connected(partner))
                ),
                None,
            )
            if bonded is not None:
                logger.warning(
                    "Ligand %s is covalently bonded to a non-terminal residue %s; "
                    "it cannot be placed in the g4composer export",
                    ligand.full_name,
                    bonded.full_name,
                )
                continue
        if bonded is not None:
            if bonded is first_by_chain[ligand.chain]:
                # ligand is the 5'-side neighbor of the first residue
                index = block_start[ligand.chain]
                insertions.append(
                    (
                        index,
                        (ligand.chain, ligand.one_letter_name, ligand, False),
                    )
                )
                block_start[ligand.chain] += 1
            else:
                index = block_end[ligand.chain]
                insertions.append(
                    (
                        index,
                        (ligand.chain, ligand.one_letter_name, ligand, False),
                    )
                )
                block_end[ligand.chain] += 1
            continue

        # base-paired ligand: place it at the sequence end opposite its partner
        pair = next(
            (
                base_pair
                for base_pair in canonical
                if (base_pair.nt1_3d is ligand and base_pair.nt2_3d in span_residues)
                or (base_pair.nt2_3d is ligand and base_pair.nt1_3d in span_residues)
            ),
            None,
        )
        if pair is None:
            continue
        partner = pair.nt2_3d if pair.nt1_3d is ligand else pair.nt1_3d
        if partner is last_by_chain[ligand.chain]:
            index = block_start[ligand.chain]
        elif partner is first_by_chain[ligand.chain]:
            index = block_end[ligand.chain]
        else:
            logger.warning(
                "Ligand %s base-pairs with non-terminal residue %s; "
                "it cannot be placed in the g4composer export",
                ligand.full_name,
                partner.full_name,
            )
            continue
        insertions.append(
            (index, (ligand.chain, ligand.one_letter_name, ligand, False))
        )
        if index == block_start[ligand.chain]:
            block_start[ligand.chain] += 1
        else:
            block_end[ligand.chain] += 1
    return insertions


def extend_entries_with_terminal_ligands(analysis: Analysis, quadruplex, entries):
    insertions = terminal_ligand_letters(analysis, quadruplex, entries)
    if not insertions:
        return entries
    result = list(entries)
    for index, insertion in sorted(insertions, key=lambda item: item[0], reverse=True):
        result.insert(index, insertion)
    return result


def unplaced_ligands(analysis: Analysis, quadruplex: Quadruplex) -> list[str]:
    """Return descriptions of nucleotide ligands that neither belong to the
    exported sequence nor were placed at the chain termini.

    Such a ligand (e.g. 4BW in 6k3y) hangs around the quadruplex chain without
    being representable in the g4composer input - it may be filling a vacant
    tetrad position that no annotator can see, so the entry cannot be trusted.
    """
    chains = export_chains(quadruplex)
    exported = {
        residue
        for _, _, residue, _ in export_entries(analysis, quadruplex)
        if residue is not None
    }
    names = []
    for residue in analysis.structure3d.residues:
        if (
            residue.is_nucleotide
            and residue.chain in chains
            and residue not in analysis.polymer_residues
            and residue not in exported
            and all(residue not in tetrad.nucleotides for tetrad in quadruplex.tetrads)
        ):
            names.append(residue.full_name)
    return names


def molecule_by_chain(entries) -> dict[str, Molecule]:
    molecules: dict[str, Molecule] = {}
    for chain, _, residue, _ in entries:
        if residue is not None and chain not in molecules:
            molecules[chain] = residue.molecule_type
    return molecules


def export_sequence(
    entries: Sequence[tuple[str, str, Residue3D | None, bool]],
    molecules: dict[str, Molecule],
) -> str:
    letters = []
    for chain, entry_letter, residue, _ in entries:
        is_rna = (
            residue.molecule_type if residue else molecules[chain]
        ) == Molecule.RNA
        letter = entry_letter.upper() if is_rna else entry_letter.lower()
        if not is_rna and letter == "u":
            # DNA residues detected by atom heuristics (e.g. thymidine
            # ligands) may carry the RNA one-letter code from rnapolis
            letter = "t"
        if letter.upper() not in "ACGTU":
            letter = "n"
        letters.append(letter)
    return "".join(letters)


def export_structure(
    analysis: Analysis,
    quadruplex: Quadruplex,
    entries: Sequence[tuple[str, str, Residue3D | None, bool]],
) -> str:
    """Overlay tetrad markers on canonical-base-pair brackets.

    g4composer uses one structure line, so we keep canonical base-pair context in
    the exported chain span and then overwrite tetrad residues with ``^`` to make
    quadruplex participation explicit.     Canonical pairs involving a tetrad
    residue are dropped entirely (on both sides) so the dot-bracket stays
    balanced and no stray brackets remain around the quadruplex core. Entries
    without coordinates (missing residues, numbering gaps) are rendered as
    ``.`` and never paired.
    """
    modelled = [residue for _, _, residue, _ in entries if residue is not None]
    tetrad_residues = {nt for tetrad in quadruplex.tetrads for nt in tetrad.nucleotides}
    brackets = canonical_dot_bracket(analysis, modelled, tetrad_residues)

    structure = []
    modelled_position = 0
    for _, _, residue, unknown in entries:
        if unknown:
            # missing / numbering-gap slots: no coordinates, rendered inert
            structure.append(".")
            continue
        char = brackets[modelled_position]
        if residue in tetrad_residues:
            char = "^"
        structure.append(char)
        modelled_position += 1
    return "".join(structure)


def canonical_dot_bracket(
    analysis: Analysis,
    residues: Sequence[Residue3D],
    tetrad_residues: set[Residue3D] | None = None,
) -> str:
    """Return a canonical-base-pair-only dot-bracket string for the export span.

    Pairs with at least one tetrad nucleotide are skipped so that the
    quadruplex core is rendered purely with ``^`` markers and the remaining
    context brackets stay balanced. A residue found in several canonical pairs
    (e.g. a bifurcated terminal pairing in 6ftu) keeps only its longest-range
    pair - the stem pairing - so every residue is bracketed at most once.
    """
    residue_to_index = {residue: i + 1 for i, residue in enumerate(residues)}
    entries = [
        Entry(i + 1, residue.one_letter_name, 0) for i, residue in enumerate(residues)
    ]

    candidates: list[tuple[int, int]] = []
    for base_pair in analysis.canonical():
        nt1 = base_pair.nt1_3d
        nt2 = base_pair.nt2_3d
        if nt1 not in residue_to_index or nt2 not in residue_to_index:
            continue
        if tetrad_residues and (nt1 in tetrad_residues or nt2 in tetrad_residues):
            continue
        i = residue_to_index[nt1]
        j = residue_to_index[nt2]
        if i == j:
            continue
        candidates.append((i, j) if i < j else (j, i))

    assigned: dict[int, int] = {}
    for i, j in sorted(
        set(candidates), key=lambda pair: pair[1] - pair[0], reverse=True
    ):
        if i in assigned or j in assigned:
            continue
        assigned[i] = j
        assigned[j] = i

    for i, j in assigned.items():
        if i < j:
            entries[i - 1].pair = j
            entries[j - 1].pair = i

    return BpSeq(entries).fcfs.structure


def export_chi(
    entries: Sequence[tuple[str, str, Residue3D | None, bool]],
) -> str:
    return "".join(
        "S" if residue is not None and residue.chi_class == GlycosidicBond.syn else "."
        for _, _, residue, _ in entries
    )


def export_sugar(
    analysis: Analysis,
    entries: Sequence[tuple[str, str, Residue3D | None, bool]],
) -> str:
    symbols = []
    for _, _, residue, _ in entries:
        if residue is None:
            symbols.append(".")
            continue
        sugar = analysis.sugar_puckers.get(residue)
        if sugar == SugarPucker.NORTH:
            symbols.append("N")
        elif sugar == SugarPucker.SOUTH:
            symbols.append("S")
        else:
            symbols.append(".")
    return "".join(symbols)


def export_orient(
    quadruplex: Quadruplex,
    build_order: Sequence[Tetrad],
    tetrad_to_label: dict[Tetrad, str],
) -> str:
    """Map tetrad polarities onto g4composer orientation signs.

    ElTetrado already exposes tetrad polarity as ``clockwise`` or
    ``anticlockwise``. g4composer uses the same concept but encodes it as ``+``
    and ``-`` respectively, emitted in g4composer build order rather than stack
    order.
    """
    polarity_by_tetrad = dict(
        zip(quadruplex.tetrads, quadruplex.tetrad_polarities, strict=True)
    )
    values = []
    for tetrad in build_order:
        polarity = polarity_by_tetrad.get(tetrad)
        if polarity is None:
            raise ValueError(
                f"g4composer export requires tetrad polarity for {tetrad_to_label[tetrad]}"
            )
        sign = "+" if polarity.value == "clockwise" else "-"
        values.append(f"{tetrad_to_label[tetrad]}{sign}")
    return ";".join(values)


def export_path(
    quadruplex: Quadruplex, mode: G4ComposerMode = G4ComposerMode.CLASSIC
) -> str:
    """Convert ElTetrado path numbering into g4composer path numbering.

    In ``CLASSIC`` mode both formats share tetrad letters and preserve ``A1`` as
    the anchor, but g4composer numbers the remaining columns clockwise whereas
    ElTetrado tracts are numbered anticlockwise.

    In ``CONSISTENT`` mode nucleotides within each tetrad are numbered 1-4 by
    the 5'->3' order in which the chain visits them (the first tetrad
    nucleotide visited by the chain is always number 1), so the clockwise
    remap does not apply.
    """
    if mode is G4ComposerMode.CONSISTENT:
        return ";".join(consistent_path(quadruplex))
    return ";".join(remap_path_entry_to_clockwise(entry) for entry in quadruplex.path)


def consistent_path(quadruplex: Quadruplex) -> list[str]:
    """Return path entries numbered per tetrad by the chain-visit order and
    the g4-axis geometry.

    The first guanine of a tetrad visited by the chain becomes number 1 and
    every other guanine of that tetrad follows in the direction that is
    clockwise when viewed along the g4 axis (the export axis, oriented so
    that the first build-order rise is positive). The numbering is therefore
    geometric, independent of both the orient and the classic tract ids.
    """
    tetrad_labels = quadruplex.tetrad_labels_by_5p_order()
    axis, centroids = quadruplex_axis_and_centroids(quadruplex)
    if axis is None:
        raise ValueError("g4composer consistent export requires the quadruplex axis")

    labels: dict[Residue3D, str] = {}
    for tetrad, tetrad_label in zip(quadruplex.tetrads, tetrad_labels, strict=True):
        visited = sorted(tetrad.nucleotides, key=lambda nt: quadruplex.global_index[nt])
        first = visited[0]
        center = centroids[tetrad]
        angles = {}
        for nt in tetrad.nucleotides:
            point = residue_c1prime_coordinates(nt)
            if point is None:
                raise ValueError(
                    f"g4composer consistent export requires C1' for {nt.full_name}"
                )
            vector = point - center
            vector = vector - numpy.dot(vector, axis) * axis
            reference = numpy.array([1.0, 0.0, 0.0])
            if abs(float(numpy.dot(reference, axis))) > 0.9:
                reference = numpy.array([0.0, 1.0, 0.0])
            reference = reference - numpy.dot(reference, axis) * axis
            angles[nt] = math.atan2(
                float(numpy.dot(numpy.cross(reference, vector), axis)),
                float(numpy.dot(reference, vector)),
            )
        # the axis-clockwise cycle: guanines by decreasing angle, rotated to
        # start at the first guanine visited by the chain
        cycle = sorted(tetrad.nucleotides, key=lambda nt: angles[nt], reverse=True)
        start = next(position for position, nt in enumerate(cycle) if nt is first)
        numbering = cycle[start:] + cycle[:start]
        for rank, nt in enumerate(numbering, start=1):
            labels[nt] = f"{tetrad_label}{rank}"
    return [
        labels[nt]
        for nt in sorted(labels.keys(), key=lambda nt: quadruplex.global_index[nt])
    ]


def remap_path_entry_to_clockwise(entry: str) -> str:
    label = "".join(ch for ch in entry if ch.isalpha())
    column = "".join(ch for ch in entry if ch.isdigit())
    if column == "1":
        return entry
    clockwise_column = {
        "2": "4",
        "3": "3",
        "4": "2",
    }.get(column)
    if clockwise_column is None:
        raise ValueError(f"Unsupported g4composer path column: {entry}")
    return f"{label}{clockwise_column}"


def export_rise(quadruplex: Quadruplex) -> str:
    return ";".join(
        format_number(rise) for rise, _ in build_order_intervals(quadruplex)
    )


def export_twist(
    quadruplex: Quadruplex, mode: G4ComposerMode = G4ComposerMode.CLASSIC
) -> str:
    if mode is G4ComposerMode.CONSISTENT:
        return ";".join(
            format_number(twist) for twist in consistent_twist_intervals(quadruplex)
        )
    return ";".join(
        format_number(twist) for _, twist in build_order_intervals(quadruplex)
    )


def consistent_twist_intervals(quadruplex: Quadruplex) -> list[float]:
    """Compute the twist for each consecutive pair in build order.

    Twist is the torsion angle C1'(1) - center(tetrad 1) - center(tetrad 2) -
    C1'(1), where ``1`` denotes the nucleotide numbered 1 in ``CONSISTENT``
    path mode (the first nucleotide of the tetrad visited by the chain),
    multiplied by the sign of the interval's rise. The sign correction
    expresses the rotation around the fixed, rise-positive quadruplex axis
    (the frame g4composer builds in): for a downward step (negative rise) the
    axis-frame rotation is the opposite of the travel-direction torsion.
    """
    build_order = tetrads_in_build_order(quadruplex)
    if len(build_order) < 2:
        return []

    axis, centroids = quadruplex_axis_and_centroids(quadruplex)
    rises: list[float] | None = None
    if axis is not None:
        rises = [
            float(numpy.dot(centroids[tetrad2] - centroids[tetrad1], axis))
            for tetrad1, tetrad2 in itertools.pairwise(build_order)
        ]
    nucleotides_by_tetrad = consistent_nucleotides_by_tetrad(quadruplex)

    intervals = []
    for interval, (tetrad1, tetrad2) in enumerate(itertools.pairwise(build_order)):
        nt1 = nucleotides_by_tetrad[tetrad1][0]
        nt2 = nucleotides_by_tetrad[tetrad2][0]
        c1 = residue_c1prime_coordinates(nt1)
        c2 = residue_c1prime_coordinates(nt2)
        if c1 is None or c2 is None:
            intervals.append(math.nan)
            continue
        twist = math.degrees(
            calculate_torsion_angle_coords(c1, tetrad1.center(), tetrad2.center(), c2)
        )
        if (
            rises is not None
            and not math.isnan(rises[interval])
            and rises[interval] != 0
        ):
            twist *= numpy.sign(rises[interval])
        intervals.append(twist)
    return intervals


def consistent_nucleotides_by_tetrad(
    quadruplex: Quadruplex,
) -> dict[Tetrad, list[Residue3D]]:
    """Map each tetrad to its nucleotides in 5'->3' chain visitation order."""
    return {
        tetrad: sorted(tetrad.nucleotides, key=lambda nt: quadruplex.global_index[nt])
        for tetrad in quadruplex.tetrads
    }


def tetrads_in_build_order(quadruplex: Quadruplex) -> list[Tetrad]:
    """Infer g4composer build order from the exported path.

    Tetrad labels are assigned by 5'->3' tetrad order, not by stack order. The
    g4composer build order is the order in which those labels first appear in the
    path, e.g. ``A1;B1;B4;A4;C4;...`` implies build order ``A, B, C``.
    """
    label_to_tetrad = {
        label: tetrad
        for tetrad, label in zip(
            quadruplex.tetrads, quadruplex.tetrad_labels_by_5p_order(), strict=True
        )
    }
    build_labels = []
    for path_entry in quadruplex.path:
        label = "".join(ch for ch in path_entry if ch.isalpha())
        if label not in build_labels:
            build_labels.append(label)
    return [label_to_tetrad[label] for label in build_labels]


def build_order_intervals(quadruplex: Quadruplex) -> list[tuple[float, float]]:
    """Compute (rise, twist) for each consecutive pair in g4composer build order.

    A single, quadruplex-wide axis (see ``quadruplex_axis_and_centroids``) is
    used for every interval so that signs are comparable across the whole
    build order, even when it reverses or bridges ElTetrado's internal
    (discovery-order) stack of adjacent tetrad pairs.
    """
    build_order = tetrads_in_build_order(quadruplex)
    if len(build_order) < 2:
        return []

    axis, centroids = quadruplex_axis_and_centroids(quadruplex)
    if axis is None:
        return [(math.nan, math.nan) for _ in build_order[1:]]

    tetrad_row_index = {
        tetrad: index for index, tetrad in enumerate(quadruplex.tetrads)
    }

    intervals = []
    for tetrad1, tetrad2 in itertools.pairwise(build_order):
        rise = float(numpy.dot(centroids[tetrad2] - centroids[tetrad1], axis))
        twist = tract_matched_twist(
            quadruplex, tetrad1, tetrad2, tetrad_row_index, axis
        )
        intervals.append((rise, twist))

    return intervals


def quadruplex_axis_and_centroids(
    quadruplex: Quadruplex,
) -> tuple[numpy.ndarray | None, dict[Tetrad, numpy.ndarray]]:
    """Return one consistently oriented axis for the whole quadruplex, plus
    each tetrad's nucleobase centroid.

    g4composer ``rise``/``twist`` are exported between consecutive tetrads in
    build order, which may differ from ElTetrado's internal (discovery-order)
    stack of adjacent tetrad pairs. Deriving a fresh axis per adjacent pair (as
    ``TetradPair`` geometry does) only yields a *locally* meaningful sign;
    chaining or reversing such per-pair values does not represent a single,
    comparable direction across the whole build order. Instead we build one
    shared axis here: each tetrad's own plane normal (arbitrarily signed by
    SVD) is chain-aligned to its discovery-order neighbor, the results are
    averaged, and the final axis is oriented so that the rise of the first
    build-order interval is positive (the description direction g4composer
    users expect).
    """
    tetrads = quadruplex.tetrads
    centroids: dict[Tetrad, numpy.ndarray] = {}
    normals: list[numpy.ndarray] = []

    for tetrad in tetrads:
        coords = collect_nucleobase_atoms(tetrad.nucleotides)
        if len(coords) < 3:
            return None, centroids
        centroid, normal = get_plane_vectors(numpy.array(coords))
        centroids[tetrad] = centroid
        normals.append(normal)

    if len(normals) < 2:
        return None, centroids

    consistent_normals = [normals[0]]
    for normal in normals[1:]:
        aligned = -normal if numpy.dot(normal, consistent_normals[-1]) < 0 else normal
        consistent_normals.append(aligned)

    axis = numpy.sum(consistent_normals, axis=0)
    norm = numpy.linalg.norm(axis)
    if norm < 1.0e-6:
        return None, centroids
    axis = axis / norm

    # Orient the axis so that the first build-order rise is positive - the
    # direction g4composer users expect the quadruplex to be described in.
    build_order = tetrads_in_build_order(quadruplex)
    if len(build_order) >= 2:
        first_rise_direction = centroids[build_order[1]] - centroids[build_order[0]]
        if numpy.dot(axis, first_rise_direction) < 0:
            axis = -axis

    return axis, centroids


def tract_matched_twist(
    quadruplex: Quadruplex,
    tetrad1: Tetrad,
    tetrad2: Tetrad,
    tetrad_row_index: dict[Tetrad, int],
    axis: numpy.ndarray,
) -> float:
    """Best-fit rotation between two build-order tetrads' C1' atoms.

    Tracts already provide a stable, quadruplex-wide nucleotide-to-nucleotide
    correspondence (one nucleotide per tetrad per tract), so this works even
    when the two tetrads are not adjacent in ElTetrado's internal stack order.
    """
    index1 = tetrad_row_index[tetrad1]
    index2 = tetrad_row_index[tetrad2]

    points1 = []
    points2 = []
    for tract in quadruplex.tracts:
        c1 = residue_c1prime_coordinates(tract.nucleotides[index1])
        c2 = residue_c1prime_coordinates(tract.nucleotides[index2])
        if c1 is None or c2 is None:
            return math.nan
        points1.append(c1)
        points2.append(c2)

    return calculate_best_fit_rotation_around_axis(
        numpy.array(points1), numpy.array(points2), axis
    )


def format_number(value: float) -> str:
    if math.isnan(value):
        return "."
    rounded = round(value, 1)
    if math.isclose(rounded, round(rounded), abs_tol=1.0e-9):
        return str(round(rounded))
    return f"{rounded:.1f}".rstrip("0").rstrip(".")


def write_g4composer(
    analysis: Analysis,
    input_path: str,
    output_path: str,
    mode: G4ComposerMode = G4ComposerMode.CLASSIC,
    allow_abasic_nucleotides: bool = False,
    allow_non_g_tetrads: bool = False,
) -> None:
    """Write one g4composer input file per involved chain.

    The ``output_path`` acts as a name template: ``/dir/4ni9-assembly1.inp``
    produces ``4ni9-assembly1-A.inp``, ``4ni9-assembly1-B.inp`` ... (a second
    distinct entry for the same chain becomes ``4ni9-assembly1-Aq2.inp``).
    Quadruplexes that raise during export are skipped; quadruplexes with a
    tetrad completed by a ligand nucleotide are skipped with a warning, as
    g4composer cannot model them. Quadruplexes containing a mixed-base tetrad
    (e.g. G·C·G·C) are skipped as well unless ``allow_non_g_tetrads`` is set -
    a partial view built of the remaining G-tetrads would misrepresent the
    structure. Entries containing unresolvable nucleotide letters (N/X) are
    skipped unless ``allow_abasic_nucleotides`` is set, in which case they are
    exported with ``n`` in the sequence. If no valid entry remains, an empty
    file is written at ``output_path``.
    """
    directory = os.path.dirname(output_path)
    root, _ = os.path.splitext(os.path.basename(output_path))

    quadruplexes = eligible_quadruplexes(analysis)
    per_chain: dict[str, list[str]] = {}
    for quad in quadruplexes:
        chains = export_chains(quad)
        if not allow_non_g_tetrads and any(
            len({nt.one_letter_name.upper() for nt in tetrad.nucleotides}) > 1
            for tetrad in quad.tetrads
        ):
            logger.warning(
                "Quadruplex contains a mixed-base tetrad; no g4composer input "
                "is generated for it unless --allow-non-g-tetrads is set"
            )
            continue
        if any(
            nt not in analysis.polymer_residues
            for tetrad in quad.tetrads
            for nt in tetrad.nucleotides
        ):
            logger.warning(
                "Quadruplex contains a tetrad completed by a ligand nucleotide; "
                "g4composer input cannot be generated for it"
            )
            continue
        try:
            entry = generate_g4composer_entry(
                analysis, quad, input_name(input_path), mode
            )
        except ValueError:
            continue
        unplaced = unplaced_ligands(analysis, quad)
        if unplaced:
            logger.warning(
                "Skipping g4composer entry for %s: unplaced ligand nucleotide(s) "
                "%s may fill a quadruplex position that g4composer cannot model",
                entry.name,
                ", ".join(unplaced),
            )
            continue
        if not allow_abasic_nucleotides:
            unresolved = sorted(
                {letter for letter in entry.sequence.upper() if letter not in "ACGTU"}
            )
            if unresolved:
                logger.warning(
                    "Skipping g4composer entry for %s: unresolvable nucleotide letter(s) %s "
                    "(use --allow-abasic-nucleotides to export them as n)",
                    entry.name,
                    ", ".join(unresolved),
                )
                continue
        blocks = per_chain.setdefault(chains[0], [])
        serialized = entry.serialize()
        if serialized in blocks:
            continue
        blocks.append(serialized)

    if not per_chain:
        with open(output_path, "w") as handle:
            handle.write("")
        return
    for chain, blocks in sorted(per_chain.items()):
        for number, block in enumerate(blocks, start=1):
            suffix = "" if number == 1 else f"q{number}"
            path = os.path.join(directory, f"{root}-{chain}{suffix}.inp")
            with open(path, "w") as handle:
                handle.write(block)
