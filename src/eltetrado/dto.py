import math
import re
from collections import Counter
from dataclasses import dataclass

from rnapolis.common import GlycosidicBond
from rnapolis.tertiary import Residue3D

from eltetrado.analysis import (
    ONZ,
    Analysis,
    LoopType,
    Quadruplex,
    Tetrad,
    TetradPair,
)
from eltetrado.model import Ion


@dataclass
class MetalDTO:
    symbol: str
    count: int


@dataclass
class NucleotideDTO:
    index: int
    chain: str
    number: int
    icode: str | None
    molecule: str
    fullName: str
    shortName: str
    chi: float
    glycosidicBond: str | None
    sugarPucker: str | None


@dataclass
class BasePairDTO:
    nt1: str
    nt2: str
    lw: str
    inTetrad: bool
    canonical: bool


@dataclass
class IonOutsideDTO:
    nt: str
    ion: str


@dataclass
class TetradDTO:
    id: str
    nt1: str
    nt2: str
    nt3: str
    nt4: str
    onz: str
    gbaClassification: str | None
    planarityDeviation: float
    ionsChannel: list[str]
    ionsOutside: list[IonOutsideDTO]


@dataclass
class TetradPairDTO:
    tetrad1: str
    tetrad2: str
    direction: str
    rise: float
    twist: float


@dataclass
class LoopDTO:
    type: str | None
    nucleotides: list[str]
    missing: list[str]


@dataclass
class SnapbackDTO:
    orientation: str
    nucleotides: list[str]
    missing: list[str]


@dataclass
class BulgeDTO:
    nucleotides: list[str]
    missing: list[str]


@dataclass
class LoopClassificationDTO:
    classification: str
    loopProgression: str


@dataclass
class QuadruplexDTO:
    tetrads: list[TetradDTO]
    onzm: str | None
    handedness: str | None
    tetradPolarities: list[str | None]
    strandPolarities: list[list[str | None]]
    loopClassification: LoopClassificationDTO | None
    gbaClassification: list[str]
    tracts: list[list[str]]
    path: list[str]
    bulges: list[BulgeDTO]
    loops: list[LoopDTO]
    snapbacks: list[SnapbackDTO]
    tags: list[str]


@dataclass
class HelixDTO:
    quadruplexes: list[QuadruplexDTO]
    tetradPairs: list[TetradPairDTO]


@dataclass
class TwoLineDotBracketDTO:
    sequence: str
    line1: str
    line2: str


@dataclass
class QuadruplexDotBracketDTO:
    sequence: str
    structure: str
    chi: str
    sugar: str
    loop: str


@dataclass
class ResultDTO:
    metals: list[MetalDTO]
    nucleotides: list[NucleotideDTO]
    basePairs: list[BasePairDTO]
    helices: list[HelixDTO]
    dotBracket: TwoLineDotBracketDTO
    quadruplexDotBracket: QuadruplexDotBracketDTO


def convert_metals(analysis: Analysis) -> list[MetalDTO]:
    counter = Counter(atom.name.title() for atom in analysis.ions)
    return [MetalDTO(Ion[k.title()].value, v) for k, v in counter.items()]


def convert_nucleotides(analysis: Analysis) -> list[NucleotideDTO]:
    return [
        NucleotideDTO(
            analysis.global_index[nt],
            nt.chain,
            nt.number,
            nt.icode,
            nt.molecule_type.value,
            nt.full_name,
            nt.one_letter_name,
            math.degrees(nt.chi),
            nt.chi_class.value if nt.chi_class else None,
            (
                analysis.sugar_puckers[nt].value
                if analysis.sugar_puckers.get(nt) is not None
                else None
            ),
        )
        for nt in analysis.structure3d.residues
        if nt.is_nucleotide
    ]


def convert_base_pairs(analysis: Analysis) -> list[BasePairDTO]:
    return [
        BasePairDTO(
            bp.nt1.full_name,
            bp.nt2.full_name,
            bp.lw.value,
            analysis.is_basepair_in_tetrad(bp),
            bp.is_canonical,
        )
        for bp in analysis.mapping.base_pairs
        if analysis.global_index[bp.nt1_3d] < analysis.global_index[bp.nt2_3d]
    ]


def convert_tetrads(quadruplex: Quadruplex) -> list[TetradDTO]:
    def id_(tetrad: Tetrad) -> str:
        return f"{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}"

    def ions_channel(tetrad: Tetrad) -> list[str]:
        return [Ion[atom.name.title()].value for atom in tetrad.ions_channel]

    def ions_outside(tetrad: Tetrad) -> list[IonOutsideDTO]:
        return [
            IonOutsideDTO(nt.full_name, Ion[atom.name.title()].value)
            for nt, atoms in tetrad.ions_outside.items()
            for atom in atoms
        ]

    return [
        TetradDTO(
            id_(tetrad),
            tetrad.nt1.full_name,
            tetrad.nt2.full_name,
            tetrad.nt3.full_name,
            tetrad.nt4.full_name,
            tetrad.onz.value,
            tetrad.gba_class.value if tetrad.gba_class is not None else None,
            float(tetrad.planarity_deviation["rmsd"]),
            ions_channel(tetrad),
            ions_outside(tetrad),
        )
        for tetrad in quadruplex.tetrads
    ]


def convert_tetrad_pairs(tetrad_pairs: list[TetradPair]) -> list[TetradPairDTO]:
    def id_(tetrad: Tetrad) -> str:
        return f"{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}"

    return [
        TetradPairDTO(
            id_(tp.tetrad1),
            id_(tp.tetrad2),
            tp.direction.value,
            float(tp.rise),
            float(tp.twist),
        )
        for tp in tetrad_pairs
    ]


def convert_quadruplexes(quadruplexes: list[Quadruplex]) -> list[QuadruplexDTO]:
    def nts_(nts: list[Residue3D]) -> list[str]:
        return [nt.full_name for nt in nts]

    return [
        QuadruplexDTO(
            convert_tetrads(q),
            q.onzm.value if q.onzm else None,
            q.handedness.value if q.handedness else None,
            [p.value if p else None for p in q.tetrad_polarities],
            [
                [p.value if p is not None else None for p in tract_polarities]
                for tract_polarities in q.strand_polarities
            ],
            (
                LoopClassificationDTO(
                    q.loop_class.value, q.loop_class.loop_progression()
                )
                if q.loop_class
                else None
            ),
            [g.value for g in q.gba_classes],
            [
                nts_(q.tracts[0].nucleotides),
                nts_(q.tracts[1].nucleotides),
                nts_(q.tracts[2].nucleotides),
                nts_(q.tracts[3].nucleotides),
            ],
            q.path,
            [
                BulgeDTO(nts_(bulge.nucleotides), list(bulge.missing))
                for bulge in q.bulges
            ],
            [
                LoopDTO(
                    loop.loop_type.value if loop.loop_type is not None else None,
                    nts_(loop.nucleotides),
                    list(loop.missing),
                )
                for loop in q.loops
            ],
            [
                SnapbackDTO(s.orientation.value, nts_(s.nucleotides), list(s.missing))
                for s in q.snapbacks
            ],
            list(q.tags),
        )
        for q in quadruplexes
    ]


def convert_helices(analysis: Analysis) -> list[HelixDTO]:
    helices = []
    for h in analysis.helices:
        # take only quadruplexes with at least two tetrads
        quadruplexes = [q for q in h.quadruplexes if len(q.tetrads) >= 2]
        # gather tetrads from selected quadruplexes
        tetrads = {t for q in quadruplexes for t in q.tetrads}
        # gather pairs among all selected tetrads
        tetrad_pairs = [
            tp
            for tp in h.tetrad_pairs
            if tp.tetrad1 in tetrads and tp.tetrad2 in tetrads
        ]

        if quadruplexes and tetrad_pairs:
            helices.append(
                HelixDTO(
                    convert_quadruplexes(quadruplexes),
                    convert_tetrad_pairs(tetrad_pairs),
                )
            )
    return helices


def convert_dot_bracket(analysis: Analysis) -> TwoLineDotBracketDTO:
    return TwoLineDotBracketDTO(analysis.sequence, analysis.line1, analysis.line2)


def convert_quadruplex_dot_bracket(analysis: Analysis) -> QuadruplexDotBracketDTO:
    letters = "QRSTUVWXYZ" + "".join(reversed("ABCDEFGHIJKLMNOP"))

    if len(analysis.tetrads) > len(letters):
        raise ValueError(f"Cannot represent more than {len(letters)} tetrads")

    lines = analysis.mapping.dot_bracket.splitlines()
    sequences = []
    structures = []

    i = 0
    while i < len(lines):
        sequences.append(lines[i + 1])
        structures.append(lines[i + 2])
        i += 3

    # Position of every modelled residue within the dashed, gap-aware
    # strings: the mapping's BpSeq covers residues and gap placeholders
    # alike, and strand lengths account for both.
    index_to_residue = analysis.mapping.bpseq_index_to_residue_map
    strands = analysis.mapping.strands_sequences
    mapping = {}
    gap_positions = set()
    bpseq_offset = 0
    dashed_offset = 0

    for (chain, strand_sequence), strand_structure in zip(
        strands, structures, strict=True
    ):
        assert chain and len(strand_sequence) == len(strand_structure)
        for j in range(len(strand_sequence)):
            residue = index_to_residue.get(bpseq_offset + j + 1)
            if residue is not None:
                mapping[residue] = dashed_offset + j
            else:
                gap_positions.add(dashed_offset + j)
        bpseq_offset += len(strand_sequence)
        dashed_offset += len(strand_sequence) + 1  # the '-' separator

    sequence = "-".join(sequences)
    structure = list("-".join(structures))
    assert len(sequence) == len(structure)
    assert len(sequence) == dashed_offset - 1

    for i, tetrad in enumerate(analysis.tetrads):
        nt1, nt2, nt3, nt4 = sorted(tetrad.nucleotides)
        if tetrad.onz == ONZ.O_PLUS:
            structure[mapping[nt1]] = letters[i]
            structure[mapping[nt2]] = letters[i].lower()
            structure[mapping[nt3]] = letters[i]
            structure[mapping[nt4]] = letters[i].lower()
        elif tetrad.onz == ONZ.O_MINUS:
            structure[mapping[nt1]] = letters[i].lower()
            structure[mapping[nt2]] = letters[i]
            structure[mapping[nt3]] = letters[i].lower()
            structure[mapping[nt4]] = letters[i]
        elif tetrad.onz == ONZ.N_PLUS:
            structure[mapping[nt1]] = letters[i]
            structure[mapping[nt2]] = letters[i].lower()
            structure[mapping[nt3]] = letters[i].lower()
            structure[mapping[nt4]] = letters[i]
        elif tetrad.onz == ONZ.N_MINUS:
            structure[mapping[nt1]] = letters[i].lower()
            structure[mapping[nt2]] = letters[i]
            structure[mapping[nt3]] = letters[i]
            structure[mapping[nt4]] = letters[i].lower()
        elif tetrad.onz == ONZ.Z_PLUS:
            structure[mapping[nt1]] = letters[i]
            structure[mapping[nt2]] = letters[i]
            structure[mapping[nt3]] = letters[i].lower()
            structure[mapping[nt4]] = letters[i].lower()
        elif tetrad.onz == ONZ.Z_MINUS:
            structure[mapping[nt1]] = letters[i].lower()
            structure[mapping[nt2]] = letters[i].lower()
            structure[mapping[nt3]] = letters[i]
            structure[mapping[nt4]] = letters[i]
        else:
            raise ValueError(f"Unexpected ONZ value: {tetrad.onz}")

    chi_line = list(re.sub(r"[^-]", ".", sequence))
    sugar_line = list(re.sub(r"[^-]", ".", sequence))
    # residues without coordinates are marked '?' in all annotation lines
    for line in (structure, chi_line, sugar_line):
        for position in gap_positions:
            line[position] = "?"

    for nt in analysis.structure3d.residues:
        if nt.is_nucleotide and nt in analysis.global_index:
            chi_line[mapping[nt]] = (
                "a"
                if nt.chi_class == GlycosidicBond.anti
                else "s"
                if nt.chi_class == GlycosidicBond.syn
                else "?"
            )
            sugar_line[mapping[nt]] = (
                "N"
                if analysis.sugar_puckers.get(nt) is not None
                and analysis.sugar_puckers[nt].value == "North"
                else "S"
                if analysis.sugar_puckers.get(nt) is not None
                and analysis.sugar_puckers[nt].value == "South"
                else "?"
            )

    loop_line = list(re.sub(r"[^-]", ".", sequence))
    for position in gap_positions:
        loop_line[position] = "?"

    for helix in analysis.helices:
        for quadruplex in helix.quadruplexes:
            for loop in quadruplex.loops:
                for nt in loop.nucleotides:
                    loop_line[mapping[nt]] = (
                        "d"
                        if loop.loop_type == LoopType.diagonal
                        else "P"
                        if loop.loop_type == LoopType.propeller_plus
                        else "p"
                        if loop.loop_type == LoopType.propeller_minus
                        else "L"
                        if loop.loop_type == LoopType.lateral_plus
                        else "l"
                        if loop.loop_type == LoopType.lateral_minus
                        else "?"
                    )
            # mark bulges
            for bulge in quadruplex.bulges:
                for nt in bulge.nucleotides:
                    loop_line[mapping[nt]] = "b"
            # mark snapbacks
            for snapback in quadruplex.snapbacks:
                for nt in snapback.nucleotides:
                    loop_line[mapping[nt]] = "s"

    return QuadruplexDotBracketDTO(
        sequence,
        "".join(structure),
        "".join(chi_line),
        "".join(sugar_line),
        "".join(loop_line),
    )


def generate_dto(analysis: Analysis):
    metals = convert_metals(analysis)
    nucleotides = convert_nucleotides(analysis)
    base_pairs = convert_base_pairs(analysis)
    helices = convert_helices(analysis)
    dot_bracket = convert_dot_bracket(analysis)
    quadruplex_dot_bracket = convert_quadruplex_dot_bracket(analysis)
    return ResultDTO(
        metals, nucleotides, base_pairs, helices, dot_bracket, quadruplex_dot_bracket
    )
