import math
import re
from collections import Counter
from dataclasses import dataclass
from typing import List, Optional

from rnapolis.common import GlycosidicBond

from eltetrado.analysis import ONZ, Analysis, LoopType, Quadruplex, TetradPair
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
    icode: Optional[str]
    molecule: str
    fullName: str
    shortName: str
    chi: float
    glycosidicBond: Optional[str]


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
    gbaClassification: Optional[str]
    planarityDeviation: float
    ionsChannel: List[str]
    ionsOutside: List[IonOutsideDTO]


@dataclass
class TetradPairDTO:
    tetrad1: str
    tetrad2: str
    direction: str
    rise: float
    twist: float


@dataclass
class LoopDTO:
    type: Optional[str]
    nucleotides: List[str]


@dataclass
class LoopClassificationDTO:
    classification: str
    loopProgression: str


@dataclass
class QuadruplexDTO:
    tetrads: List[TetradDTO]
    onzm: Optional[str]
    loopClassification: Optional[LoopClassificationDTO]
    gbaClassification: List[str]
    tracts: List[List[str]]
    bulges: List[str]
    loops: List[LoopDTO]


@dataclass
class HelixDTO:
    quadruplexes: List[QuadruplexDTO]
    tetradPairs: List[TetradPairDTO]


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
    loop: str


@dataclass
class ResultDTO:
    metals: List[MetalDTO]
    nucleotides: List[NucleotideDTO]
    basePairs: List[BasePairDTO]
    helices: List[HelixDTO]
    dotBracket: TwoLineDotBracketDTO
    quadruplexDotBracket: QuadruplexDotBracketDTO


def convert_metals(analysis: Analysis) -> List[MetalDTO]:
    counter = Counter(map(lambda atom: atom.name.title(), analysis.ions))
    return [MetalDTO(Ion[k.title()].value, v) for k, v in counter.items()]


def convert_nucleotides(analysis: Analysis) -> List[NucleotideDTO]:
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
        )
        for nt in analysis.structure3d.residues
        if nt.is_nucleotide
    ]


def convert_base_pairs(analysis: Analysis) -> List[BasePairDTO]:
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


def convert_tetrads(quadruplex: Quadruplex) -> List[TetradDTO]:
    id_ = (
        lambda tetrad: f"{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}"
    )
    ions_channel = lambda tetrad: [
        Ion[atom.name.title()].value for atom in tetrad.ions_channel
    ]
    ions_outside = lambda tetrad: [
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
            float(tetrad.planarity_deviation),
            ions_channel(tetrad),
            ions_outside(tetrad),
        )
        for tetrad in quadruplex.tetrads
    ]


def convert_tetrad_pairs(tetrad_pairs: List[TetradPair]) -> List[TetradPairDTO]:
    id_ = (
        lambda tetrad: f"{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}"
    )
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


def convert_quadruplexes(quadruplexes: List[Quadruplex]) -> List[QuadruplexDTO]:
    nts_ = lambda nts: [nt.full_name for nt in nts]
    return [
        QuadruplexDTO(
            convert_tetrads(q),
            q.onzm.value if q.onzm else None,
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
            [nt.full_name for nt in q.bulges],
            [
                LoopDTO(
                    l.loop_type.value if l.loop_type is not None else None,
                    nts_(l.nucleotides),
                )
                for l in q.loops
            ],
        )
        for q in quadruplexes
    ]


def convert_helices(analysis: Analysis) -> List[HelixDTO]:
    helices = []
    for h in analysis.helices:
        # take only quadruplexes with at least two tetrads
        quadruplexes = [q for q in h.quadruplexes if len(q.tetrads) >= 2]
        # gather tetrads from selected quadruplexes
        tetrads = set([t for q in quadruplexes for t in q.tetrads])
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

    mapping = {}
    i = 0
    current_chain = None

    for residue in analysis.structure3d.residues:
        if residue.is_nucleotide:
            if current_chain is None:
                current_chain = residue.chain
            if residue.chain != current_chain:
                current_chain = residue.chain
                i += 1
            mapping[residue] = i
            i += 1

    sequence = "-".join(sequences)
    structure = list("-".join(structures))
    assert len(sequence) == len(structure)
    assert len(sequence) == i

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

    for nt in analysis.structure3d.residues:
        if nt.is_nucleotide and nt in analysis.global_index:
            chi_line[mapping[nt]] = (
                "a"
                if nt.chi_class == GlycosidicBond.anti
                else "s"
                if nt.chi_class == GlycosidicBond.syn
                else "?"
            )

    loop_line = list(re.sub(r"[^-]", ".", sequence))

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
            for nt in quadruplex.bulges:
                loop_line[mapping[nt]] = "b"

    return QuadruplexDotBracketDTO(
        sequence, "".join(structure), "".join(chi_line), "".join(loop_line)
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
