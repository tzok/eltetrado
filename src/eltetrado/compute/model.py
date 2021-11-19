from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import List

import numpy

METALS = {
    x.casefold() for x in [
        'Ag', 'Au', 'Ba', 'Ca', 'Co', 'Cs', 'Cu', 'Eu', 'Fe', 'Ir', 'K', 'Li', 'Mg', 'Mn', 'Na', 'Ni', 'Os', 'Pb', 'Pt',
        'Ru', 'Sr', 'Tl', 'V', 'Zn'
    ]
}


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
    other = '...'


class Stericity(Enum):
    cis = 'cis'
    trans = 'trans'


class Edge(Enum):
    WC = 'Watson-Crick'
    H = 'Hoogsteen'
    S = 'Sugar'


class ONZ(Enum):
    O_PLUS = 'O+'
    O_MINUS = 'O-'
    N_PLUS = 'N+'
    N_MINUS = 'N-'
    Z_PLUS = 'Z+'
    Z_MINUS = 'Z-'


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


class LoopType(Enum):
    diagonal = 'diagonal'
    propeller_plus = 'propeller+'
    propeller_minus = 'propeller-'
    lateral_plus = 'lateral+'
    lateral_minus = 'lateral-'


class Direction(Enum):
    parallel = 'parallel'
    antiparallel = 'antiparallel'
    hybrid = 'hybrid'


@dataclass
class Metal:
    symbol: Ion
    count: int


@dataclass
class Nucleotide:
    index: int
    model: int
    chain: str
    number: int
    icode: str
    molecule: Molecule
    fullName: str
    shortName: str
    chi: float
    glycosidicBond: GlycosidicBond


@dataclass
class BasePair:
    nt1: str
    nt2: str
    stericity: Stericity
    edge5: Edge
    edge3: Edge


@dataclass
class IonOutside:
    nt: str
    ion: Ion


@dataclass
class Tetrad:
    id: str
    nt1: str
    nt2: str
    nt3: str
    nt4: str
    onz: ONZ
    gbaClassification: GbaTetradClassification
    planarityDeviation: float
    ionsChannel: List[Ion]
    ionsOutside: List[IonOutside]


@dataclass
class Loop:
    type: LoopType
    nucleotides: List[str]


@dataclass
class Quadruplex:
    tetrads: List[Tetrad]
    onzm: ONZM
    loopClassification: LoopClassification
    gbaClassification: List[GbaQuadruplexClassification]
    tracts: List[List[str]]
    loops: List[Loop]


@dataclass
class TetradPair:
    tetrad1: str
    tetrad2: str
    direction: Direction
    rise: float
    twist: float


@dataclass
class Helix:
    quadruplexes: List[Quadruplex]
    tetradPairs: List[TetradPair]


@dataclass
class TwoLineDotBracket:
    sequence: str
    line1: str
    line2: str


@dataclass
class Result:
    metals: List[Metal]
    nucleotides: List[Nucleotide]
    basePairs: List[BasePair]
    helices: List[Helix]
    dotBracket: TwoLineDotBracket


@dataclass
class Atom3D:
    atom_name: str
    residue_name: str
    chain_identifier: str
    residue_number: int
    insertion_code: str
    x: float
    y: float
    z: float

    def coordinates(self):
        return numpy.array([self.x, self.y, self.z])


@dataclass
class Residue3D:
    atoms: List[Atom3D]
    residue_name: str
    chain_identifier: str
    residue_number: int
    insertion_code: str


@dataclass
class Structure3D:
    residues: List[Residue3D]

    def find_metal_ions(self):
        atoms = []
        used = set()
        for residue in self.residues:
            for atom in residue.atoms:
                if atom.atom_name.casefold() in METALS and tuple(atom.coordinates()) not in used:
                    atoms.append(atom)
                    used.add(tuple(atom.coordinates()))
        return atoms


def convert_metals(analysis):
    counter = Counter(map(lambda atom: atom.atom_name.title(), analysis.metal_ions))
    return [Metal(Ion(k.title()), v) for k, v in counter.items()]


def convert_nucleotides(analysis):
    return [
        Nucleotide(nt.index, nt.model, nt.chain, nt.number, nt.icode,
                   Molecule(nt.molecule), nt.full_name, nt.short_name, float(nt.chi if nt.chi else 'nan'),
                   GlycosidicBond(nt.glycosidic_bond)) for nt in analysis.nucleotides.values()
    ]


def convert_base_pairs(analysis):
    lw2stericity = lambda lw: Stericity.cis if lw[0] == 'c' else Stericity.trans
    lw2edge = lambda lw: Edge.WC if lw == 'W' else Edge.H if lw == 'H' else Edge.S
    return [
        BasePair(bp.pair[0].full_name, bp.pair[1].full_name, lw2stericity(bp.lw), lw2edge(bp.lw[1]), lw2edge(bp.lw[2]))
        for bp in analysis.pairs.values()
    ]


def convert_ions_channel(tetrad):
    return [Ion(atom.atom_name.title()) for atom in tetrad.ions_channel]


def convert_ions_outside(tetrad):
    return [
        IonOutside(nt.full_name, Ion(ion.atom_name.title()))
        for nt, ions_list in tetrad.ions_outside.items()
        for ion in ions_list
    ]


def convert_tetrads(quadruplex):
    return [
        Tetrad(repr(t), t.nucleotides[0].full_name, t.nucleotides[1].full_name,
               t.nucleotides[2].full_name, t.nucleotides[3].full_name, ONZ(t.get_classification()),
               GbaTetradClassification(t.gba_classification()), float(t.planarity_deviation), convert_ions_channel(t),
               convert_ions_outside(t)) for t in quadruplex.tetrads
    ]


def convert_tracts(quadruplex):
    return [[nt.full_name for nt in tract.nucleotides] for tract in quadruplex.tracts]


def convert_loops(quadruplex):
    return [Loop(LoopType(loop.loop_type), [nt.full_name for nt in loop.nucleotides]) for loop in quadruplex.loops]


def convert_quadruplexes(helix):
    return [
        Quadruplex(convert_tetrads(q), ONZM(f'{q.onzm_classification()}{q.direction()}{q.sign()}'),
                   LoopClassification(q.loop_classification),
                   [GbaQuadruplexClassification(gba) for gba in q.gba_classification], convert_tracts(q),
                   convert_loops(q)) for q in filter(lambda q: len(q.tetrads) > 1, helix.quadruplexes)
    ]


def convert_tetrad_pairs(helix):
    return [
        TetradPair(repr(tp.tetrad1), repr(tp.tetrad2), Direction(tp.direction), float(tp.rise), float(tp.twist))
        for tp in helix.tetrad_pairs
    ]


def convert_helices(analysis):
    return [
        Helix(convert_quadruplexes(h), convert_tetrad_pairs(h))
        for h in filter(lambda h: len(h.tetrads) > 1, analysis.helices)
    ]


def convert_dot_bracket(analysis):
    return TwoLineDotBracket(analysis.sequence, analysis.line1, analysis.line2)


def generate_dto(analysis):
    metals = convert_metals(analysis)
    nucleotides = convert_nucleotides(analysis)
    base_pairs = convert_base_pairs(analysis)
    helices = convert_helices(analysis)
    dot_bracket = convert_dot_bracket(analysis)
    return Result(metals, nucleotides, base_pairs, helices, dot_bracket)
