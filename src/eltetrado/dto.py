import math
from collections import Counter
from dataclasses import dataclass
from typing import List, Optional

from eltetrado.analysis import Analysis, Helix, Quadruplex
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
class ResultDTO:
    metals: List[MetalDTO]
    nucleotides: List[NucleotideDTO]
    basePairs: List[BasePairDTO]
    helices: List[HelixDTO]
    dotBracket: TwoLineDotBracketDTO


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
    ions_channel = lambda tetrad: [atom.name for atom in tetrad.ions_channel]
    ions_outside = lambda tetrad: [
        IonOutsideDTO(nt.full_name, atom.name)
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


def convert_tetrad_pairs(helix: Helix) -> List[TetradPairDTO]:
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
        for tp in helix.tetrad_pairs
    ]


def convert_quadruplexes(helix: Helix) -> List[QuadruplexDTO]:
    nts_ = lambda nts: [nt.full_name for nt in nts]
    return [
        QuadruplexDTO(
            convert_tetrads(q),
            q.onzm.value if q.onzm else None,
            LoopClassificationDTO(q.loop_class.value, q.loop_class.loop_progression())
            if q.loop_class
            else None,
            [g.value for g in q.gba_classes],
            [
                nts_(q.tracts[0].nucleotides),
                nts_(q.tracts[1].nucleotides),
                nts_(q.tracts[2].nucleotides),
                nts_(q.tracts[3].nucleotides),
            ],
            [
                LoopDTO(
                    l.loop_type.value if l.loop_type is not None else None,
                    nts_(l.nucleotides),
                )
                for l in q.loops
            ],
        )
        for q in helix.quadruplexes
    ]


def convert_helices(analysis: Analysis) -> List[HelixDTO]:
    """
    NOTE: there is an "if" which will prevent single-tetrad helices from serialization; this is on purpose
    """
    return [
        HelixDTO(convert_quadruplexes(h), convert_tetrad_pairs(h))
        for h in analysis.helices
        if len(h.tetrads) > 1
    ]


def convert_dot_bracket(analysis: Analysis) -> TwoLineDotBracketDTO:
    return TwoLineDotBracketDTO(analysis.sequence, analysis.line1, analysis.line2)


def generate_dto(analysis: Analysis):
    metals = convert_metals(analysis)
    nucleotides = convert_nucleotides(analysis)
    base_pairs = convert_base_pairs(analysis)
    helices = convert_helices(analysis)
    dot_bracket = convert_dot_bracket(analysis)
    return ResultDTO(metals, nucleotides, base_pairs, helices, dot_bracket)
