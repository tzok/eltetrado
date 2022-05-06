import logging
import math
import string
from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import Enum
from typing import List, Optional, Dict, Tuple

import numpy


class LeontisWesthof(Enum):
    cWW = 'cWW'
    cWH = 'cWH'
    cWS = 'cWS'
    cHW = 'cHW'
    cHH = 'cHH'
    cHS = 'cHS'
    cSW = 'cSW'
    cSH = 'cSH'
    cSS = 'cSS'
    tWW = 'tWW'
    tWH = 'tWH'
    tWS = 'tWS'
    tHW = 'tHW'
    tHH = 'tHH'
    tHS = 'tHS'
    tSW = 'tSW'
    tSH = 'tSH'
    tSS = 'tSS'


class Saenger(Enum):
    I = 'I'
    II = 'II'
    III = 'III'
    IV = 'IV'
    V = 'V'
    VI = 'VI'
    VII = 'VII'
    VIII = 'VIII'
    IX = 'IX'
    X = 'X'
    XI = 'XI'
    XII = 'XII'
    XIII = 'XIII'
    XIV = 'XIV'
    XV = 'XV'
    XVI = 'XVI'
    XVII = 'XVII'
    XVIII = 'XVIII'
    XIX = 'XIX'
    XX = 'XX'
    XXI = 'XXI'
    XXII = 'XXII'
    XXIII = 'XXIII'
    XXIV = 'XXIV'
    XXV = 'XXV'
    XXVI = 'XXVI'
    XXVII = 'XXVII'
    XXVIII = 'XXVIII'


class StackingTopology(Enum):
    upward = 'upward'
    downward = 'downward'
    inward = 'inward'
    outward = 'outward'


class BR(Enum):
    _0 = 0
    _1 = 1
    _2 = 2
    _3 = 3
    _4 = 4
    _5 = 5
    _6 = 6
    _7 = 7
    _8 = 8
    _9 = 9


class BPh(Enum):
    _0 = 0
    _1 = 1
    _2 = 2
    _3 = 3
    _4 = 4
    _5 = 5
    _6 = 6
    _7 = 7
    _8 = 8
    _9 = 9


@dataclass(frozen=True)
class ResidueLabel:
    chain: str
    number: int
    name: str


@dataclass(frozen=True)
class ResidueAuth:
    chain: str
    number: int
    icode: str
    name: str


@dataclass
class Residue:
    label: Optional[ResidueLabel]
    auth: Optional[ResidueAuth]

    def __post_init__(self):
        if isinstance(self.label, dict):
            self.label = ResidueLabel(**self.label)
        if isinstance(self.auth, dict):
            self.auth = ResidueAuth(**self.auth)


@dataclass
class Interaction:
    nt1: Residue
    nt2: Residue

    def __post_init__(self):
        if isinstance(self.nt1, dict):
            self.nt1 = Residue(**self.nt1)
        if isinstance(self.nt2, dict):
            self.nt2 = Residue(**self.nt2)


@dataclass
class BasePair(Interaction):
    lw: LeontisWesthof
    saenger: Optional[Saenger]

    def __post_init__(self):
        super(BasePair, self).__post_init__()
        if isinstance(self.lw, str):
            self.lw = LeontisWesthof[self.lw]
        if isinstance(self.saenger, str):
            self.saenger = Saenger[self.saenger]


@dataclass
class Stacking(Interaction):
    topology: Optional[StackingTopology]

    def __post_init__(self):
        super(Stacking, self).__post_init__()
        if isinstance(self.topology, str):
            self.topology = StackingTopology[self.topology]


@dataclass
class BaseRibose(Interaction):
    br: Optional[BR]

    def __post_init__(self):
        super(BaseRibose, self).__post_init__()
        if isinstance(self.br, str):
            self.br = BR[self.br]


@dataclass
class BasePhosphate(Interaction):
    bph: Optional[BPh]

    def __post_init__(self):
        super(BasePhosphate, self).__post_init__()
        if isinstance(self.bph, str):
            self.bph = BPh[self.bph]


@dataclass
class OtherInteraction(Interaction):
    def __post_init__(self):
        super(OtherInteraction, self).__post_init__()


@dataclass
class Structure2D:
    basePairs: List[BasePair]
    stackings: List[Stacking]
    baseRiboseInteractions: List[BaseRibose]
    basePhosphateInteractions: List[BasePhosphate]
    otherInteractions: List[OtherInteraction]

    def __post_init__(self):
        self.basePairs = [BasePair(**x) for x in self.basePairs if isinstance(x, dict)]
        self.stackings = [Stacking(**x) for x in self.stackings if isinstance(x, dict)]
        self.baseRiboseInteractions = [BaseRibose(**x) for x in self.baseRiboseInteractions if isinstance(x, dict)]
        self.basePhosphateInteractions = \
            [BasePhosphate(**x) for x in self.basePhosphateInteractions if isinstance(x, dict)]
        self.otherInteractions = [OtherInteraction(**x) for x in self.otherInteractions if isinstance(x, dict)]


class Ion(Enum):
    Ag = 'Ag'
    Au = 'Au'
    Ba = 'Ba'
    Ca = 'Ca'
    Co = 'Co'
    Cs = 'Cs'
    Cu = 'Cu'
    Eu = 'Eu'
    Fe = 'Fe'
    Ir = 'Ir'
    K = 'K'
    Li = 'Li'
    Mg = 'Mg'
    Mn = 'Mn'
    Na = 'Na'
    Ni = 'Ni'
    Os = 'Os'
    Pb = 'Pb'
    Pt = 'Pt'
    Ru = 'Ru'
    Sr = 'Sr'
    Tl = 'Tl'
    V = 'V'
    Zn = 'Zn'


class Molecule(Enum):
    DNA = 'DNA'
    RNA = 'RNA'
    Other = 'Other'


class GlycosidicBond(Enum):
    anti = 'anti'
    syn = 'syn'


class ONZ(Enum):
    O_PLUS = 'O+'
    O_MINUS = 'O-'
    N_PLUS = 'N+'
    N_MINUS = 'N-'
    Z_PLUS = 'Z+'
    Z_MINUS = 'Z-'

    def score(self):
        if self == ONZ.O_PLUS:
            return 0
        elif self == ONZ.O_MINUS:
            return 1
        elif self == ONZ.N_PLUS:
            return 2
        elif self == ONZ.N_MINUS:
            return 3
        elif self == ONZ.Z_PLUS:
            return 4
        elif self == ONZ.Z_MINUS:
            return 5
        return 6


class ONZM(Enum):
    OP_PLUS = 'Op+'
    OP_MINUS = 'Op-'
    OP_STAR = 'Op*'
    OA_PLUS = 'Oa+'
    OA_MINUS = 'Oa-'
    OA_STAR = 'Oa*'
    OH_PLUS = 'Oh+'
    OH_MINUS = 'Oh-'
    OH_STAR = 'Oh*'
    NP_PLUS = 'Np+'
    NP_MINUS = 'Np-'
    NP_STAR = 'Np*'
    NA_PLUS = 'Na+'
    NA_MINUS = 'Na-'
    NA_STAR = 'Na*'
    NH_PLUS = 'Nh+'
    NH_MINUS = 'Nh-'
    NH_STAR = 'Nh*'
    ZP_PLUS = 'Zp+'
    ZP_MINUS = 'Zp-'
    ZP_STAR = 'Zp*'
    ZA_PLUS = 'Za+'
    ZA_MINUS = 'Za-'
    ZA_STAR = 'Za*'
    ZH_PLUS = 'Zh+'
    ZH_MINUS = 'Zh-'
    ZH_STAR = 'Zh*'
    MP_PLUS = 'Mp+'
    MP_MINUS = 'Mp-'
    MP_STAR = 'Mp*'
    MA_PLUS = 'Ma+'
    MA_MINUS = 'Ma-'
    MA_STAR = 'Ma*'
    MH_PLUS = 'Mh+'
    MH_MINUS = 'Mh-'
    MH_STAR = 'Mh*'

    @staticmethod
    def from_value(value: str):
        for candidate in ONZM:
            if candidate.value == value:
                return candidate
        raise RuntimeError(f'Failed to match {value} to an ONZM class')


class GbaTetradClassification(Enum):
    Ia = 'Ia'
    IIa = 'IIa'
    IIIa = 'IIIa'
    IVa = 'IVa'
    Va = 'Va'
    VIa = 'VIa'
    VIIa = 'VIIa'
    VIIIa = 'VIIIa'
    Ib = 'Ib'
    IIb = 'IIb'
    IIIb = 'IIIb'
    IVb = 'IVb'
    Vb = 'Vb'
    VIb = 'VIb'
    VIIb = 'VIIb'
    VIIIb = 'VIIIb'


class GbaQuadruplexClassification(Enum):
    I = 'I'
    II = 'II'
    III = 'III'
    IV = 'IV'
    V = 'V'
    VI = 'VI'
    VII = 'VII'
    VIII = 'VIII'


class LoopClassification(Enum):
    _1a = '1a'
    _2a = '2a'
    _3a = '3a'
    _4a = '4a'
    _5a = '5a'
    _6a = '6a'
    _7a = '7a'
    _8a = '8a'
    _9a = '9a'
    _10a = '10a'
    _11a = '11a'
    _12a = '12a'
    _13a = '13a'
    _1b = '1b'
    _2b = '2b'
    _3b = '3b'
    _4b = '4b'
    _5b = '5b'
    _6b = '6b'
    _7b = '7b'
    _8b = '8b'
    _9b = '9b'
    _10b = '10b'
    _11b = '11b'
    _12b = '12b'
    _13b = '13b'
    invalid = 'n/a'

    def loop_progression(self) -> str:
        if self == LoopClassification._1a:
            return '-(ppp)'
        elif self == LoopClassification._1b:
            return '+(ppp)'
        elif self == LoopClassification._2a:
            return '-(ppl)'
        elif self == LoopClassification._2b:
            return '+(ppl)'
        elif self == LoopClassification._3a:
            return '-(plp)'
        elif self == LoopClassification._3b:
            return '+(plp)'
        elif self == LoopClassification._4a:
            return '-(lpp)'
        elif self == LoopClassification._4b:
            return '+(lpp)'
        elif self == LoopClassification._5a:
            return '(-pd+p)'
        elif self == LoopClassification._5b:
            return '(+pd-p)'
        elif self == LoopClassification._6a:
            return '-(lll)'
        elif self == LoopClassification._6b:
            return '+(lll)'
        elif self == LoopClassification._7a:
            return '-(llp)'
        elif self == LoopClassification._7b:
            return '+(llp)'
        elif self == LoopClassification._8a:
            return '-(lpl)'
        elif self == LoopClassification._8b:
            return '+(lpl)'
        elif self == LoopClassification._9a:
            return '-(pll)'
        elif self == LoopClassification._9b:
            return '+(pll)'
        elif self == LoopClassification._10a:
            return '(-pd+l)'
        elif self == LoopClassification._10b:
            return '(+pd-l)'
        elif self == LoopClassification._11a:
            return '(-ld+l)'
        elif self == LoopClassification._11b:
            return '(+ld-l)'
        elif self == LoopClassification._12a:
            return '(d-pd)'
        elif self == LoopClassification._12b:
            return '(d+pd)'
        elif self == LoopClassification._13a:
            return '(-ld+p)'
        elif self == LoopClassification._13b:
            return '(+ld-p)'
        raise RuntimeError(f'Failed to get string representation of {self}')

    @staticmethod
    def from_value(value: str):
        for candidate in LoopClassification:
            if candidate.value == value:
                return candidate
        raise RuntimeError(f'Failed to match {value} to a LoopClassification class')


class LoopType(Enum):
    diagonal = 'diagonal'
    propeller_plus = 'propeller+'
    propeller_minus = 'propeller-'
    lateral_plus = 'lateral+'
    lateral_minus = 'lateral-'

    @staticmethod
    def from_value(value: str):
        for loop_type in LoopType:
            if loop_type.value == value:
                return loop_type
        raise RuntimeError(f'Failed to match {value} to a LoopType class')


class Direction(Enum):
    parallel = 'parallel'
    antiparallel = 'antiparallel'
    hybrid = 'hybrid'


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
    gbaClassification: str
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
    type: str
    nucleotides: List[str]


@dataclass
class LoopClassificationDTO:
    classification: str
    loop_progression: str


@dataclass
class QuadruplexDTO:
    tetrads: List[TetradDTO]
    onzm: str
    loopClassification: LoopClassificationDTO
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


@dataclass(frozen=True)
class Atom3D:
    label: Optional[ResidueLabel]
    auth: Optional[ResidueAuth]
    model: int
    atomName: str
    x: float
    y: float
    z: float

    def coordinates(self) -> numpy.ndarray:
        return numpy.array([self.x, self.y, self.z])


@dataclass(order=True)
class Residue3D:
    index: int
    name: str
    one_letter_name: str
    model: int
    label: Optional[ResidueLabel]
    auth: Optional[ResidueAuth]
    atoms: Tuple[Atom3D]

    # Dict representing expected name of atom involved in glycosidic bond
    outermost_atoms = {'A': 'N9', 'G': 'N9', 'C': 'N1', 'U': 'N1', 'T': 'N1'}
    # Dist representing expected name of atom closest to the tetrad center
    innermost_atoms = {'A': 'N6', 'G': 'O6', 'C': 'N4', 'U': 'O4', 'T': 'O4'}

    def __hash__(self):
        return hash((self.name, self.model, self.label, self.auth, self.atoms))

    @property
    def chain(self) -> str:
        if self.auth is not None:
            return self.auth.chain
        return self.label.chain

    @property
    def number(self) -> int:
        if self.auth is not None:
            return self.auth.number
        return self.label.number

    @property
    def icode(self) -> Optional[str]:
        if self.auth is not None:
            return self.auth.icode if self.auth.icode not in (' ', '?') else None
        return None

    @property
    def molecule_type(self) -> Molecule:
        if self.name.upper() in ('A', 'C', 'G', 'U'):
            return Molecule.RNA
        if self.name.upper() in ('DA', 'DC', 'DG', 'DT'):
            return Molecule.DNA
        return Molecule.Other

    @property
    def full_name(self) -> str:
        if self.auth:
            builder = f'{self.auth.chain}.{self.auth.name}'
            if self.auth.name[-1] in string.digits:
                builder += '/'
            builder += f'{self.auth.number}'
            if self.auth.icode not in (' ', '?'):
                builder += f'^{self.auth.icode}'
            return builder
        else:
            builder = f'{self.label.chain}.{self.label.name}'
            if self.label.name[-1] in string.digits:
                builder += '/'
            builder += f'{self.label.number}'
            return builder

    @property
    def chi(self) -> float:
        if self.one_letter_name.upper() in ('A', 'G'):
            return self.__chi_purine()
        elif self.one_letter_name.upper() in ('C', 'U', 'T'):
            return self.__chi_pyrimidine()
        # if unknown, try purine first, then pyrimidine
        torsion = self.__chi_purine()
        if math.isnan(torsion):
            return self.__chi_pyrimidine()
        return torsion

    # TODO: the ranges could be modified to match Saenger
    @property
    def chi_class(self) -> Optional[GlycosidicBond]:
        if math.isnan(self.chi):
            return None
        if -math.pi / 2 < self.chi < math.pi / 2:
            return GlycosidicBond.syn
        return GlycosidicBond.anti

    def find_atom(self, atom_name: str) -> Optional[Atom3D]:
        for atom in self.atoms:
            if atom.atomName == atom_name:
                return atom
        return None

    @property
    def outermost_atom(self) -> Atom3D:
        return next(filter(None, self.__outer_generator()))

    @property
    def innermost_atom(self) -> Atom3D:
        return next(filter(None, self.__inner_generator()))

    @property
    def is_nucleotide(self) -> bool:
        return len(self.atoms) > 1 and any([atom for atom in self.atoms if atom.atomName == "C1'"])

    def __chi_purine(self) -> float:
        atoms = [self.find_atom("O4'"), self.find_atom("C1'"), self.find_atom("N9"), self.find_atom("C4")]
        if all(atoms):
            return Residue3D.__torsion_angle(atoms)
        return math.nan

    def __chi_pyrimidine(self) -> float:
        atoms = [self.find_atom("O4'"), self.find_atom("C1'"), self.find_atom("N1"), self.find_atom("C2")]
        if all(atoms):
            return Residue3D.__torsion_angle(atoms)
        return math.nan

    @staticmethod
    def __torsion_angle(atoms: List[Atom3D]) -> float:
        assert len(atoms) == 4
        v1 = atoms[1].coordinates() - atoms[0].coordinates()
        v2 = atoms[2].coordinates() - atoms[1].coordinates()
        v3 = atoms[3].coordinates() - atoms[2].coordinates()
        t1 = numpy.cross(v1, v2)
        t2 = numpy.cross(v2, v3)
        t3 = v1 * numpy.linalg.norm(v2)
        return math.atan2(numpy.dot(t2, t3), numpy.dot(t1, t2))

    def __outer_generator(self):
        # try to find expected atom name
        upper = self.one_letter_name.upper()
        if upper in self.outermost_atoms:
            yield self.find_atom(self.outermost_atoms[upper])

        # try to get generic name for purine/pyrimidine
        yield self.find_atom('N9')
        yield self.find_atom('N1')

        # try to find at least C1' next to nucleobase
        yield self.find_atom("C1'")

        # get any atom
        if self.atoms:
            yield self.atoms[0]

        # last resort, create pseudoatom at (0, 0, 0)
        logging.error(
            f'Failed to determine the outermost atom for nucleotide {self}, so an arbitrary atom will be used')
        yield Atom3D(self.label, self.auth, self.model, 'UNK', 0.0, 0.0, 0.0)

    def __inner_generator(self):
        # try to find expected atom name
        upper = self.one_letter_name.upper()
        if upper in self.innermost_atoms:
            yield self.find_atom(self.innermost_atoms[upper])

        # try to get generic name for purine/pyrimidine
        yield self.find_atom('C6')
        yield self.find_atom('C4')

        # try to find any atom at position 4 or 6 for purine/pyrimidine respectively
        yield self.find_atom('O6')
        yield self.find_atom('N6')
        yield self.find_atom('S6')
        yield self.find_atom('O4')
        yield self.find_atom('N4')
        yield self.find_atom('S4')

        # get any atom
        if self.atoms:
            yield self.atoms[0]

        # last resort, create pseudoatom at (0, 0, 0)
        logging.error(
            f'Failed to determine the innermost atom for nucleotide {self}, so an arbitrary atom will be used')
        yield Atom3D(self.label, self.auth, self.model, 'UNK', 0.0, 0.0, 0.0)


@dataclass(frozen=True)
class BasePair3D:
    nt1: Residue3D
    nt2: Residue3D
    lw: LeontisWesthof

    score_table = {
        LeontisWesthof.cWW: 1, LeontisWesthof.tWW: 2, LeontisWesthof.cWH: 3, LeontisWesthof.tWH: 4,
        LeontisWesthof.cWS: 5, LeontisWesthof.tWS: 6, LeontisWesthof.cHW: 7, LeontisWesthof.tHW: 8,
        LeontisWesthof.cHH: 9, LeontisWesthof.tHH: 10, LeontisWesthof.cHS: 11, LeontisWesthof.tHS: 12,
        LeontisWesthof.cSW: 13, LeontisWesthof.tSW: 14, LeontisWesthof.cSH: 15, LeontisWesthof.tSH: 16,
        LeontisWesthof.cSS: 17, LeontisWesthof.tSS: 18
    }

    def conflicts_with(self, other) -> bool:
        xi, yi = sorted([self.nt1.index, self.nt2.index])
        xj, yj = sorted([other.nt1.index, other.nt2.index])
        return xi < xj < yi < yj or xj < xi < yj < yi

    def reverse(self):
        lw = f'{self.lw.name[0]}{self.lw.name[2]}{self.lw.name[1]}'
        return BasePair3D(self.nt2, self.nt1, LeontisWesthof[lw])

    def score(self) -> int:
        return self.score_table.get(self.lw, 20)

    def is_canonical(self) -> bool:
        nts = ''.join(sorted([self.nt1.one_letter_name.upper(), self.nt2.one_letter_name.upper()]))
        return self.lw == 'cWW' and (nts == 'AU' or nts == 'CG' or nts == 'GU')


@dataclass
class Stacking3D:
    nt1: Residue3D
    nt2: Residue3D


@dataclass
class Structure3D:
    residues: List[Residue3D]

    def find_residue(self, label: Optional[ResidueLabel], auth: Optional[ResidueAuth]) -> Optional[Residue3D]:
        for residue in self.residues:
            if label is not None and label == residue.label:
                return residue
            if auth is not None and auth == residue.auth:
                return residue
        return None

    def base_pairs(self, structure2d: Structure2D) -> List[BasePair3D]:
        result = []
        for base_pair in structure2d.basePairs:
            nt1 = self.find_residue(base_pair.nt1.label, base_pair.nt1.auth)
            nt2 = self.find_residue(base_pair.nt2.label, base_pair.nt2.auth)
            result.append(BasePair3D(nt1, nt2, base_pair.lw))
        return result

    def base_pair_graph(self, structure2d: Structure2D, strict: bool = False) -> Dict[Residue3D, List[Residue3D]]:
        graph = defaultdict(list)
        for pair in self.base_pairs(structure2d):
            if strict and pair.lw not in (LeontisWesthof.cWH, LeontisWesthof.cHW):
                continue
            graph[pair.nt1].append(pair.nt2)
        return graph

    def base_pair_dict(self, structure2d: Structure2D, strict: bool = False) \
            -> Dict[Tuple[Residue3D, Residue3D], BasePair3D]:
        result = {}
        for pair in self.base_pairs(structure2d):
            if strict and pair.lw not in (LeontisWesthof.cWH, LeontisWesthof.cHW):
                continue
            result[(pair.nt1, pair.nt2)] = pair
        return result

    def stackings(self, structure2d: Structure2D) -> List[Stacking3D]:
        result = []
        for stacking in structure2d.stackings:
            nt1 = self.find_residue(stacking.nt1.label, stacking.nt1.auth)
            nt2 = self.find_residue(stacking.nt2.label, stacking.nt2.auth)
            result.append(Stacking3D(nt1, nt2))
        return result

    def stacking_graph(self, structure2d: Structure2D) -> Dict[Residue3D, List[Residue3D]]:
        graph = defaultdict(list)
        for pair in self.stackings(structure2d):
            graph[pair.nt1].append(pair.nt2)
        return graph

    def stacking_dict(self, structure2d: Structure2D) -> Dict[Tuple[Residue3D, Residue3D], Stacking3D]:
        result = {}
        for pair in self.stackings(structure2d):
            result[(pair.nt1, pair.nt2)] = pair
        return result


def convert_metals(analysis) -> List[MetalDTO]:
    counter = Counter(map(lambda atom: atom.atomName.title(), analysis.ions))
    return [MetalDTO(Ion[k.title()].value, v) for k, v in counter.items()]


def convert_nucleotides(analysis) -> List[NucleotideDTO]:
    return [
        NucleotideDTO(nt.index, nt.chain, nt.number, nt.icode, nt.molecule_type.value, nt.full_name,
                      nt.one_letter_name, math.degrees(nt.chi), nt.chi_class.value if nt.chi_class else None)
        for nt in analysis.structure3d.residues
    ]


def convert_base_pairs(analysis) -> List[BasePairDTO]:
    return [BasePairDTO(bp.nt1.full_name, bp.nt2.full_name, bp.lw.value) for bp in analysis.base_pairs]


def convert_tetrads(quadruplex) -> List[TetradDTO]:
    id_ = lambda tetrad: f'{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}'
    ions_channel = lambda tetrad: [atom.atomName for atom in tetrad.ions_channel]
    ions_outside = lambda tetrad: [IonOutsideDTO(nt.full_name, atom.atomName)
                                   for nt, atoms in tetrad.ions_outside.items() for atom in atoms]
    return [
        TetradDTO(id_(tetrad), tetrad.nt1.full_name, tetrad.nt2.full_name, tetrad.nt3.full_name, tetrad.nt4.full_name,
                  tetrad.onz.value, tetrad.gba_class.value, float(tetrad.planarity_deviation), ions_channel(tetrad),
                  ions_outside(tetrad))
        for tetrad in quadruplex.tetrads
    ]


def convert_tetrad_pairs(helix) -> List[TetradPairDTO]:
    id_ = lambda tetrad: f'{tetrad.nt1.full_name}-{tetrad.nt2.full_name}-{tetrad.nt3.full_name}-{tetrad.nt4.full_name}'
    return [
        TetradPairDTO(id_(tp.tetrad1), id_(tp.tetrad2), tp.direction.value, float(tp.rise), float(tp.twist))
        for tp in helix.tetrad_pairs
    ]


def convert_quadruplexes(helix) -> List[QuadruplexDTO]:
    nts_ = lambda nts: [nt.full_name for nt in nts]
    return [
        QuadruplexDTO(convert_tetrads(q), q.onzm.value if q.onzm else None,
                      LoopClassificationDTO(q.loop_class.value,
                                            q.loop_class.loop_progression()) if q.loop_class else None,
                      [g.value for g in q.gba_classes],
                      [nts_(q.tracts[0].nucleotides), nts_(q.tracts[1].nucleotides),
                       nts_(q.tracts[2].nucleotides), nts_(q.tracts[3].nucleotides)],
                      [LoopDTO(l.loop_type.value, nts_(l.nucleotides)) for l in q.loops])
        for q in helix.quadruplexes
    ]


def convert_helices(analysis) -> List[HelixDTO]:
    return [
        HelixDTO(convert_quadruplexes(h), convert_tetrad_pairs(h)) for h in analysis.helices
    ]


def convert_dot_bracket(analysis) -> TwoLineDotBracketDTO:
    return TwoLineDotBracketDTO(analysis.sequence, analysis.line1, analysis.line2)


def generate_dto(analysis):
    metals = convert_metals(analysis)
    nucleotides = convert_nucleotides(analysis)
    base_pairs = convert_base_pairs(analysis)
    helices = convert_helices(analysis)
    dot_bracket = convert_dot_bracket(analysis)
    return ResultDTO(metals, nucleotides, base_pairs, helices, dot_bracket)
