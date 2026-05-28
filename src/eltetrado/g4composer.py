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
"""

import math
import os
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy
from rnapolis.common import Entry, GlycosidicBond, Molecule, BpSeq
from rnapolis.tertiary import Residue3D

from eltetrado.analysis import (
    Analysis,
    Quadruplex,
    Tetrad,
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


def select_single_quadruplex(analysis: Analysis) -> Quadruplex:
    quadruplexes = eligible_quadruplexes(analysis)
    if not quadruplexes:
        raise ValueError(
            "g4composer export requires exactly one unimolecular quadruplex with at least 2 tetrads, found none"
        )
    if len(quadruplexes) > 1:
        raise ValueError(
            f"g4composer export requires exactly one unimolecular quadruplex with at least 2 tetrads, found {len(quadruplexes)}"
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
        structure=export_structure(analysis, quadruplex, residues),
        chi=export_chi(residues),
        sugar=export_sugar(analysis, residues),
        orient=export_orient(quadruplex, build_order, tetrad_to_label),
        rise=export_rise(quadruplex),
        twist=export_twist(quadruplex),
        path=export_path(quadruplex),
    )


def export_residues(analysis: Analysis, quadruplex: Quadruplex) -> List[Residue3D]:
    """Export full nucleotide chains participating in the selected quadruplex.

    g4composer expects the quadruplex together with its 5' and 3' flanking
    single-stranded context, so we export every nucleotide from the involved
    chain(s), not only residues that belong to tetrads.
    """
    chains = export_chains(quadruplex)
    return [
        residue
        for residue in sorted(
            analysis.structure3d.residues,
            key=lambda residue: analysis.global_index[residue],
        )
        if residue.is_nucleotide and residue.chain in chains
    ]


def export_chains(quadruplex: Quadruplex) -> List[str]:
    return sorted(
        {
            nt.chain
            for tetrad in quadruplex.tetrads
            for nt in tetrad.nucleotides
            if nt.is_nucleotide
        }
    )


def export_sequence(residues: Sequence[Residue3D]) -> str:
    letters = []
    for residue in residues:
        letter = residue.one_letter_name
        if residue.molecule_type == Molecule.RNA:
            letters.append(letter.upper())
        else:
            letters.append(letter.lower())
    return "".join(letters)


def export_structure(
    analysis: Analysis, quadruplex: Quadruplex, residues: Sequence[Residue3D]
) -> str:
    """Overlay tetrad markers on canonical-base-pair brackets.

    g4composer uses one structure line, so we keep canonical base-pair context in
    the exported chain span and then overwrite tetrad residues with ``^`` to make
    quadruplex participation explicit.
    """
    structure = list(canonical_dot_bracket(analysis, residues))
    tetrad_residues = {nt for tetrad in quadruplex.tetrads for nt in tetrad.nucleotides}
    for i, residue in enumerate(residues):
        if residue in tetrad_residues:
            structure[i] = "^"
    return "".join(structure)


def canonical_dot_bracket(analysis: Analysis, residues: Sequence[Residue3D]) -> str:
    """Return a canonical-base-pair-only dot-bracket string for the export span."""
    residue_to_index = {residue: i + 1 for i, residue in enumerate(residues)}
    entries = [
        Entry(i + 1, residue.one_letter_name, 0) for i, residue in enumerate(residues)
    ]
    seen = set()

    for base_pair in analysis.canonical():
        nt1 = base_pair.nt1_3d
        nt2 = base_pair.nt2_3d
        if nt1 not in residue_to_index or nt2 not in residue_to_index:
            continue
        key = tuple(sorted((residue_to_index[nt1], residue_to_index[nt2])))
        if key in seen:
            continue
        seen.add(key)
        i = residue_to_index[nt1]
        j = residue_to_index[nt2]
        entries[i - 1].pair = j
        entries[j - 1].pair = i

    return BpSeq(entries).fcfs.structure


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
    """Map tetrad polarities onto g4composer orientation signs.

    ElTetrado already exposes tetrad polarity as ``clockwise`` or
    ``anticlockwise``. g4composer uses the same concept but encodes it as ``+``
    and ``-`` respectively, emitted in g4composer build order rather than stack
    order.
    """
    polarity_by_tetrad = {
        tetrad: polarity
        for tetrad, polarity in zip(quadruplex.tetrads, quadruplex.tetrad_polarities)
    }
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


def export_path(quadruplex: Quadruplex) -> str:
    """Convert ElTetrado path numbering into g4composer path numbering.

    Both formats share tetrad letters and preserve ``A1`` as the anchor, but
    g4composer numbers the remaining columns clockwise whereas ElTetrado tracts
    are numbered anticlockwise.
    """
    return ";".join(remap_path_entry_to_clockwise(entry) for entry in quadruplex.path)


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
    intervals = build_intervals(quadruplex)
    values = [
        format_number(sum(step.signed_rise for step in interval))
        for interval in intervals
    ]
    return ";".join(values)


def export_twist(quadruplex: Quadruplex) -> str:
    intervals = build_intervals(quadruplex)
    values = [
        format_number(sum(step.signed_twist for step in interval))
        for interval in intervals
    ]
    return ";".join(values)


@dataclass(frozen=True)
class SignedStep:
    signed_rise: float
    signed_twist: float


def build_intervals(quadruplex: Quadruplex) -> List[List[SignedStep]]:
    """Build signed traversal intervals between consecutive build-order tetrads.

    g4composer ``rise`` and ``twist`` are exported between consecutive tetrads in
    build order (A->B, B->C, ...), but the underlying geometry in ElTetrado is
    defined on adjacent stacked tetrad pairs. When build order differs from stack
    order, an interval may therefore traverse one or more stack edges, possibly
    in reverse. Reverse traversal negates the adjacent-pair metric.
    """
    build_order = tetrads_in_build_order(quadruplex)
    if len(build_order) < 2:
        return []

    pair_steps = signed_pair_steps(quadruplex)
    intervals = []
    for start, end in zip(build_order, build_order[1:]):
        intervals.append(path_steps_between(start, end, pair_steps))
    return intervals


def tetrads_in_build_order(quadruplex: Quadruplex) -> List[Tetrad]:
    """Infer g4composer build order from the exported path.

    Tetrad labels are assigned by 5'->3' tetrad order, not by stack order. The
    g4composer build order is the order in which those labels first appear in the
    path, e.g. ``A1;B1;B4;A4;C4;...`` implies build order ``A, B, C``.
    """
    label_to_tetrad = {
        label: tetrad
        for tetrad, label in zip(
            quadruplex.tetrads, quadruplex.tetrad_labels_by_5p_order()
        )
    }
    build_labels = []
    for path_entry in quadruplex.path:
        label = "".join(ch for ch in path_entry if ch.isalpha())
        if label not in build_labels:
            build_labels.append(label)
    return [label_to_tetrad[label] for label in build_labels]


def signed_pair_steps(
    quadruplex: Quadruplex,
) -> Dict[Tuple[Tetrad, Tetrad], SignedStep]:
    """Cache signed adjacent-stack metrics for both traversal directions."""
    steps = {}
    for pair in quadruplex.tetrad_pairs:
        rise = signed_rise_for_pair(pair.tetrad1, pair.tetrad2)
        twist = pair.twist
        steps[(pair.tetrad1, pair.tetrad2)] = SignedStep(rise, twist)
        steps[(pair.tetrad2, pair.tetrad1)] = SignedStep(-rise, -twist)
    return steps


def path_steps_between(
    start: Tetrad,
    end: Tetrad,
    pair_steps: Dict[Tuple[Tetrad, Tetrad], SignedStep],
) -> List[SignedStep]:
    """Traverse the stack graph between two build-order tetrads.

    The returned step list is later summed to produce g4composer ``rise`` and
    ``twist`` for a single build-order interval.
    """
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
                (
                    neighbor,
                    steps + [pair_steps[(current, neighbor)]],
                    visited + [neighbor],
                )
            )

    raise ValueError("Failed to derive g4composer build traversal between tetrads")


def adjacency(edges: Iterable[Tuple[Tetrad, Tetrad]]) -> Dict[Tetrad, List[Tetrad]]:
    result: Dict[Tetrad, List[Tetrad]] = {}
    for left, right in edges:
        result.setdefault(left, []).append(right)
    return result


def signed_rise_for_pair(tetrad1: Tetrad, tetrad2: Tetrad) -> float:
    """Return the signed rise for one adjacent stacked tetrad pair."""
    coords1 = numpy.array(collect_nucleobase_atoms(tetrad1.nucleotides))
    coords2 = numpy.array(collect_nucleobase_atoms(tetrad2.nucleotides))
    if len(coords1) == 0 or len(coords2) == 0:
        return math.nan
    return calculate_signed_rise(coords1, coords2)


def format_number(value: float) -> str:
    if math.isnan(value):
        return "."
    rounded = round(value, 1)
    if math.isclose(rounded, round(rounded), abs_tol=1.0e-9):
        return str(int(round(rounded)))
    return f"{rounded:.1f}".rstrip("0").rstrip(".")


def write_g4composer(
    analysis: Analysis,
    input_path: str,
    output_path: str,
    quadruplex: Optional[Quadruplex] = None,
) -> None:
    selected = (
        quadruplex if quadruplex is not None else select_single_quadruplex(analysis)
    )
    entry = generate_g4composer_entry(analysis, selected, input_name(input_path))
    with open(output_path, "w") as handle:
        handle.write(entry.serialize())
