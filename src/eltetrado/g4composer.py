import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy
from rnapolis.common import GlycosidicBond, Molecule
from rnapolis.tertiary import Residue3D

from eltetrado.analysis import (
    Analysis,
    Quadruplex,
    Tetrad,
    calculate_quadruplex_twist_centroids,
    calculate_signed_rise,
    collect_nucleobase_atoms,
)
from eltetrado.model import SugarPucker


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


def eligible_quadruplexes(analysis: Analysis) -> List[Quadruplex]:
    return [
        quadruplex
        for helix in analysis.helices
        for quadruplex in helix.quadruplexes
        if len(quadruplex.tetrads) >= 2
    ]


def select_single_quadruplex(analysis: Analysis) -> Quadruplex:
    quadruplexes = eligible_quadruplexes(analysis)
    if not quadruplexes:
        raise ValueError(
            "g4composer export requires exactly one quadruplex with at least 2 tetrads, found none"
        )
    if len(quadruplexes) > 1:
        raise ValueError(
            f"g4composer export requires exactly one quadruplex with at least 2 tetrads, found {len(quadruplexes)}"
        )
    return quadruplexes[0]


def generate_g4composer_entry(
    analysis: Analysis, quadruplex: Quadruplex, name: str
) -> G4ComposerEntry:
    residues = export_residues(analysis, quadruplex)
    build_order = tetrads_in_build_order(quadruplex)
    tetrad_labels = quadruplex.tetrad_labels_by_5p_order()
    tetrad_to_label = {
        tetrad: tetrad_labels[index] for index, tetrad in enumerate(quadruplex.tetrads)
    }

    return G4ComposerEntry(
        name=name,
        sequence=export_sequence(residues),
        structure=export_structure(quadruplex, residues),
        chi=export_chi(residues),
        sugar=export_sugar(analysis, residues),
        orient=export_orient(quadruplex, build_order, tetrad_to_label),
        rise=export_rise(quadruplex),
        twist=export_twist(quadruplex),
        path=";".join(quadruplex.path),
    )


def export_residues(analysis: Analysis, quadruplex: Quadruplex) -> List[Residue3D]:
    chains = {
        nt.chain for tetrad in quadruplex.tetrads for nt in tetrad.nucleotides if nt.is_nucleotide
    }
    return [
        residue
        for residue in sorted(
            analysis.structure3d.residues, key=lambda residue: analysis.global_index[residue]
        )
        if residue.is_nucleotide and residue.chain in chains
    ]


def export_sequence(residues: Sequence[Residue3D]) -> str:
    letters = []
    for residue in residues:
        letter = residue.one_letter_name
        if residue.molecule_type == Molecule.RNA:
            letters.append(letter.upper())
        else:
            letters.append(letter.lower())
    return "".join(letters)


def export_structure(quadruplex: Quadruplex, residues: Sequence[Residue3D]) -> str:
    tetrad_residues = {nt for tetrad in quadruplex.tetrads for nt in tetrad.nucleotides}
    return "".join("^" if residue in tetrad_residues else "." for residue in residues)


def export_chi(residues: Sequence[Residue3D]) -> str:
    return "".join(
        "S" if residue.chi_class == GlycosidicBond.syn else "." for residue in residues
    )


def export_sugar(analysis: Analysis, residues: Sequence[Residue3D]) -> str:
    symbols = []
    for residue in residues:
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
    tetrad_to_label: Dict[Tetrad, str],
) -> str:
    centers_and_normals = {}
    for index, tetrad in enumerate(quadruplex.tetrads):
        geometry = quadruplex.tetrad_geometry(tetrad)
        if geometry is None:
            centers_and_normals[tetrad] = (None, None)
            continue
        center, normal = geometry
        reference_nt = quadruplex.tracts[0].nucleotides[index]
        normal = quadruplex.oriented_tetrad_normal(index, tetrad, normal, reference_nt)
        centers_and_normals[tetrad] = (center, normal)

    if len(build_order) > 1:
        first_center = centers_and_normals[build_order[0]][0]
        second_center = centers_and_normals[build_order[1]][0]
        stacking_direction = (
            second_center - first_center
            if first_center is not None and second_center is not None
            else None
        )
    else:
        stacking_direction = None

    values = []
    for tetrad in build_order:
        _, normal = centers_and_normals[tetrad]
        sign = "+"
        if normal is not None and stacking_direction is not None:
            sign = "+" if numpy.dot(normal, stacking_direction) >= 0 else "-"
        values.append(f"{tetrad_to_label[tetrad]}{sign}")
    return ";".join(values)


def export_rise(quadruplex: Quadruplex) -> str:
    intervals = build_intervals(quadruplex)
    values = [format_number(sum(step.signed_rise for step in interval)) for interval in intervals]
    return ";".join(values)


def export_twist(quadruplex: Quadruplex) -> str:
    intervals = build_intervals(quadruplex)
    values = [format_number(sum(step.signed_twist for step in interval)) for interval in intervals]
    return ";".join(values)


@dataclass(frozen=True)
class SignedStep:
    signed_rise: float
    signed_twist: float


def build_intervals(quadruplex: Quadruplex) -> List[List[SignedStep]]:
    build_order = tetrads_in_build_order(quadruplex)
    if len(build_order) < 2:
        return []

    pair_steps = signed_pair_steps(quadruplex)
    intervals = []
    for start, end in zip(build_order, build_order[1:]):
        intervals.append(path_steps_between(start, end, pair_steps))
    return intervals


def tetrads_in_build_order(quadruplex: Quadruplex) -> List[Tetrad]:
    label_to_tetrad = {
        label: tetrad
        for tetrad, label in zip(quadruplex.tetrads, quadruplex.tetrad_labels_by_5p_order())
    }
    build_labels = []
    for path_entry in quadruplex.path:
        label = "".join(ch for ch in path_entry if ch.isalpha())
        if label not in build_labels:
            build_labels.append(label)
    return [label_to_tetrad[label] for label in build_labels]


def signed_pair_steps(quadruplex: Quadruplex) -> Dict[Tuple[Tetrad, Tetrad], SignedStep]:
    steps = {}
    for pair in quadruplex.tetrad_pairs:
        rise = signed_rise_for_pair(pair.tetrad1, pair.tetrad2)
        twist = signed_twist_for_pair(pair.tetrad1, pair.tetrad2)
        steps[(pair.tetrad1, pair.tetrad2)] = SignedStep(rise, twist)
        steps[(pair.tetrad2, pair.tetrad1)] = SignedStep(-rise, -twist)
    return steps


def path_steps_between(
    start: Tetrad,
    end: Tetrad,
    pair_steps: Dict[Tuple[Tetrad, Tetrad], SignedStep],
) -> List[SignedStep]:
    graph = adjacency(pair_steps.keys())
    queue: List[Tuple[Tetrad, List[SignedStep], List[Tetrad]]] = [(start, [], [start])]

    while queue:
        current, steps, visited = queue.pop(0)
        if current == end:
            return steps
        for neighbor in graph.get(current, []):
            if neighbor in visited:
                continue
            queue.append(
                (neighbor, steps + [pair_steps[(current, neighbor)]], visited + [neighbor])
            )

    raise ValueError("Failed to derive g4composer build traversal between tetrads")


def adjacency(edges: Iterable[Tuple[Tetrad, Tetrad]]) -> Dict[Tetrad, List[Tetrad]]:
    result: Dict[Tetrad, List[Tetrad]] = {}
    for left, right in edges:
        result.setdefault(left, []).append(right)
    return result


def signed_rise_for_pair(tetrad1: Tetrad, tetrad2: Tetrad) -> float:
    coords1 = numpy.array(collect_nucleobase_atoms(tetrad1.nucleotides))
    coords2 = numpy.array(collect_nucleobase_atoms(tetrad2.nucleotides))
    if len(coords1) == 0 or len(coords2) == 0:
        return math.nan
    return calculate_signed_rise(coords1, coords2)


def signed_twist_for_pair(tetrad1: Tetrad, tetrad2: Tetrad) -> float:
    tetrad1_all_coords = collect_nucleobase_atoms(tetrad1.nucleotides)
    tetrad2_all_coords = collect_nucleobase_atoms(tetrad2.nucleotides)
    if not tetrad1_all_coords or not tetrad2_all_coords:
        return math.nan

    nt_list1 = [numpy.array(collect_nucleobase_atoms((nt,))) for nt in tetrad1.nucleotides]
    nt_list2 = [numpy.array(collect_nucleobase_atoms((nt,))) for nt in tetrad2.nucleotides]
    if any(len(coords) == 0 for coords in nt_list1 + nt_list2):
        return math.nan

    return float(
        calculate_quadruplex_twist_centroids(
            numpy.array(tetrad1_all_coords),
            numpy.array(tetrad2_all_coords),
            nt_list1,
            nt_list2,
        )
    )


def format_number(value: float) -> str:
    if math.isnan(value):
        return "."
    rounded = round(value, 1)
    if math.isclose(rounded, round(rounded), abs_tol=1.0e-9):
        return str(int(round(rounded)))
    return f"{rounded:.1f}".rstrip("0").rstrip(".")


def write_g4composer(
    analysis: Analysis, input_path: str, output_path: str, quadruplex: Optional[Quadruplex] = None
) -> None:
    selected = quadruplex if quadruplex is not None else select_single_quadruplex(analysis)
    entry = generate_g4composer_entry(analysis, selected, input_name(input_path))
    with open(output_path, "w") as handle:
        handle.write(entry.serialize())
