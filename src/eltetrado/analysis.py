import itertools
import logging
import math
import os
import string
import subprocess
import tempfile
from collections import defaultdict, Counter
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Tuple, Optional, Set

import numpy

from eltetrado.model import Atom3D, Structure3D, Structure2D, BasePair3D, Residue3D, GlycosidicBond, ONZ, \
    GbaTetradClassification, Ion, Direction, LoopType, ONZM, GbaQuadruplexClassification, LoopClassification

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))


@dataclass(order=True)
class Tetrad:
    @staticmethod
    def is_valid(nt1: Residue3D, nt2: Residue3D, nt3: Residue3D, nt4: Residue3D,
                 pair_dictionary: Dict[Tuple[Residue3D, Residue3D], BasePair3D]) -> bool:
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
    onz: ONZ = field(init=False)
    gba_class: Optional[GbaTetradClassification] = field(init=False)
    planarity_deviation: float = field(init=False)
    ions_channel: List[Atom3D] = field(default_factory=list)
    ions_outside: Dict[Residue3D, List[Atom3D]] = field(default_factory=dict)

    def __post_init__(self):
        self.reorder_to_match_5p_3p()
        self.planarity_deviation = self.__calculate_planarity_deviation()

    def reorder_to_match_5p_3p(self):
        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = map(lambda nt: nt.index, self.nucleotides)
        indices = sorted((ni, nj, nk, nl))
        ni, nj, nk, nl = (indices.index(x) for x in (ni, nj, nk, nl))

        nmin = min(ni, nj, nk, nl)
        if nmin == ni:
            pass
        elif nmin == nj:
            self.nt1, self.nt2, self.nt3, self.nt4 = self.nt2, self.nt3, self.nt4, self.nt1
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_23, self.pair_34, self.pair_41, self.pair_12
        elif nmin == nk:
            self.nt1, self.nt2, self.nt3, self.nt4 = self.nt3, self.nt4, self.nt1, self.nt2
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_34, self.pair_41, self.pair_12, self.pair_23
        else:
            self.nt1, self.nt2, self.nt3, self.nt4 = self.nt4, self.nt1, self.nt2, self.nt3
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_41, self.pair_12, self.pair_23, self.pair_34

        # flip order if necessary
        if self.pair_12.score() > self.pair_41.reverse().score():
            self.nt1, self.nt2, self.nt3, self.nt4 = self.nt1, self.nt4, self.nt3, self.nt2
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_41.reverse(), self.pair_34.reverse(), self.pair_23.reverse(), self.pair_12.reverse()

        # ONZ and da Silva's classification are valid in 5'-3' order
        self.onz = self.__classify_onz()
        self.gba_class = self.__classify_by_gba()

    def reorder_to_match_other_tetrad(self, order: Tuple[Residue3D, Residue3D, Residue3D, Residue3D]):
        if order == (self.nt1, self.nt2, self.nt3, self.nt4):
            pass
        elif order == (self.nt2, self.nt3, self.nt4, self.nt1):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_23, self.pair_34, self.pair_41, self.pair_12
        elif order == (self.nt3, self.nt4, self.nt1, self.nt2):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_34, self.pair_41, self.pair_12, self.pair_23
        elif order == (self.nt4, self.nt1, self.nt2, self.nt3):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_41, self.pair_12, self.pair_23, self.pair_34
        elif order == (self.nt4, self.nt3, self.nt2, self.nt1):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_34.reverse(), self.pair_23.reverse(), self.pair_12.reverse(), self.pair_41.reverse()
        elif order == (self.nt3, self.nt2, self.nt1, self.nt4):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_23.reverse(), self.pair_12.reverse(), self.pair_41.reverse(), self.pair_34.reverse()
        elif order == (self.nt2, self.nt1, self.nt4, self.nt3):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_12.reverse(), self.pair_41.reverse(), self.pair_34.reverse(), self.pair_23.reverse()
        elif order == (self.nt1, self.nt4, self.nt3, self.nt2):
            self.pair_12, self.pair_23, self.pair_34, self.pair_41 = self.pair_41.reverse(), self.pair_34.reverse(), self.pair_23.reverse(), self.pair_12.reverse()
        else:
            raise RuntimeError(f'Cannot apply order: {order}')

        self.nt1, self.nt2, self.nt3, self.nt4 = order

    def __classify_onz(self) -> ONZ:
        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = (nt.index for nt in self.nucleotides)
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

        raise RuntimeError(f'Impossible combination: {ni} {nj} {nk} {nl}')

    def __classify_by_gba(self) -> Optional[GbaTetradClassification]:
        """
        See: Webba da Silva, M. (2007). Geometric Formalism for DNA Quadruplex Folding.
        Chemistry - A European Journal, 13(35), 9738–9745. https://doi.org/10.1002/chem.200701255

        :return: Classification according to Webba da Silva or n/a
        """
        # without all nucleotides having a valid syn/anti, this classification is impossible
        if not all([nt.chi_class in (GlycosidicBond.syn, GlycosidicBond.anti) for nt in self.nucleotides]):
            return None

        # this will create a 4-letter string made of 's' for syn or 'a' for anti
        fingerprint = ''.join([nt.chi_class.value[0] for nt in self.nucleotides])

        # this dict has all classes mapped to fingerprints
        gba_classes = {
            'aass': GbaTetradClassification.Ia,
            'ssaa': GbaTetradClassification.Ib,
            'asas': GbaTetradClassification.IIa,
            'sasa': GbaTetradClassification.IIb,
            'asaa': GbaTetradClassification.IIIa,
            'sass': GbaTetradClassification.IIIb,
            'aaas': GbaTetradClassification.IVa,
            'sssa': GbaTetradClassification.IVb,
            'aasa': GbaTetradClassification.Va,
            'ssas': GbaTetradClassification.Vb,
            'assa': GbaTetradClassification.VIa,
            'saas': GbaTetradClassification.VIb,
            'asss': GbaTetradClassification.VIIa,
            'saaa': GbaTetradClassification.VIIb,
            'aaaa': GbaTetradClassification.VIIIa,
            'ssss': GbaTetradClassification.VIIIb
        }

        if fingerprint not in gba_classes:
            logging.error(f'Impossible combination of syn/anti: {[nt.chi_class for nt in self.nucleotides]}')
            return None
        return gba_classes[fingerprint]

    def __calculate_planarity_deviation(self) -> float:
        outer = [nt.outermost_atom for nt in self.nucleotides]
        inner = [nt.innermost_atom for nt in self.nucleotides]
        return numpy.linalg.norm(center_of_mass(outer) - center_of_mass(inner))

    @property
    def nucleotides(self) -> Tuple[Residue3D, Residue3D, Residue3D, Residue3D]:
        return self.nt1, self.nt2, self.nt3, self.nt4

    def __hash__(self):
        return hash(frozenset([self.nt1, self.nt2, self.nt3, self.nt4]))

    def __str__(self):
        return f'    ' \
               f'{self.nt1.full_name} {self.nt2.full_name} {self.nt3.full_name} {self.nt4.full_name} ' \
               f'{self.pair_12.lw.value} {self.pair_23.lw.value} {self.pair_34.lw.value} {self.pair_41.lw.value} ' \
               f'{self.onz.value} {self.gba_class.value} ' \
               f'planarity={round(self.planarity_deviation, 2)} ' \
               f'{self.__ions_channel_str()} ' \
               f'{self.__ions_outside_str()}\n'

    def chains(self) -> Set[str]:
        return set([nt.chain for nt in self.nucleotides])

    def is_disjoint(self, other) -> bool:
        return frozenset(self.nucleotides).isdisjoint(frozenset(other.nucleotides))

    def center(self) -> numpy.ndarray:
        return center_of_mass(self.outer_and_inner_atoms())

    def outer_and_inner_atoms(self) -> List[Atom3D]:
        return list(map(lambda residue: residue.outermost_atom, self.nucleotides)) + \
               list(map(lambda residue: residue.innermost_atom, self.nucleotides))

    def __ions_channel_str(self) -> str:
        if self.ions_channel:
            return 'ions_channel=' + ','.join([atom.atomName for atom in self.ions_channel])
        return ''

    def __ions_outside_str(self) -> str:
        if self.ions_outside:
            result = []
            for residue, ions in self.ions_outside.items():
                result.append(f'{residue.full_name}: [{",".join([ion.atomName for ion in ions])}]')
            return 'ions_outside=' + ' '.join(result)
        return ''


@dataclass
class TetradPair:
    tetrad1: Tetrad
    tetrad2: Tetrad
    stacked: Dict[Residue3D, Residue3D]
    tetrad2_nts_best_order: Tuple[Residue3D, Residue3D, Residue3D, Residue3D] = field(init=False)
    direction: Direction = field(init=False)
    rise: float = field(init=False)
    twist: float = field(init=False)

    def __post_init__(self):
        self.tetrad2_nts_best_order = (
            self.stacked[self.tetrad1.nt1], self.stacked[self.tetrad1.nt2],
            self.stacked[self.tetrad1.nt3], self.stacked[self.tetrad1.nt4]
        )
        self.direction = self.__determine_direction()
        self.rise = self.__calculate_rise()
        self.twist = self.__calculate_twist()

    def __determine_direction(self) -> Direction:
        indices1 = list(map(lambda nt: nt.index, self.tetrad1.nucleotides))
        indices2 = list(map(lambda nt: nt.index, self.tetrad2_nts_best_order))

        # count directions 5' -> 3' as +1 or -1
        counter = Counter(1 if j - i > 0 else -1 for i, j in zip(indices1, indices2))
        direction, count = counter.most_common()[0]

        if count == 4:
            # all in the same direction
            return Direction.parallel
        elif count == 2:
            # two in +, one in - direction
            return Direction.antiparallel

        return Direction.hybrid

    def __calculate_rise(self) -> float:
        t1 = self.tetrad1.outer_and_inner_atoms()
        t2 = self.tetrad2.outer_and_inner_atoms()
        return numpy.linalg.norm(center_of_mass(t1) - center_of_mass(t2))

    def __calculate_twist(self) -> float:
        nt1_1, nt1_2, _, _ = self.tetrad1.nucleotides
        nt2_1, nt2_2, _, _ = self.tetrad2_nts_best_order

        v1 = nt1_1.find_atom("C1'").coordinates() - nt1_2.find_atom("C1'").coordinates()
        v1 = v1 / numpy.linalg.norm(v1)
        v2 = nt2_1.find_atom("C1'").coordinates() - nt2_2.find_atom("C1'").coordinates()
        v2 = v2 / numpy.linalg.norm(v2)
        return math.degrees(numpy.arccos(numpy.clip(numpy.dot(v1, v2), -1.0, 1.0)))

    def __str__(self):
        return f'      direction={self.direction.value} rise={round(self.rise, 2)} twist={round(self.twist, 2)}\n'


@dataclass
class Tract:
    nucleotides: List[Residue3D]

    def __str__(self):
        return f'      {", ".join(map(lambda nt: nt.full_name, self.nucleotides))}'


@dataclass
class Loop:
    nucleotides: List[Residue3D]
    loop_type: Optional[LoopType]

    def __str__(self):
        return f'      {self.loop_type.value if self.loop_type else "n/a"} ' \
               f'{", ".join(map(lambda nt: nt.full_name, self.nucleotides))}'


@dataclass
class Quadruplex:
    tetrads: List[Tetrad]
    tetrad_pairs: List[TetradPair]
    structure3d: Structure3D
    onzm: Optional[ONZM] = field(init=False)
    gba_classes: List[GbaQuadruplexClassification] = field(init=False)
    tracts: List[Tract] = field(init=False)
    loops: List[Loop] = field(init=False)
    loop_class: Optional[LoopClassification] = field(init=False)

    def __post_init__(self):
        self.onzm = self.__classify_onzm()
        self.gba_classes = self.__classify_by_gba()
        self.tracts = self.__find_tracts()
        self.loops = self.__find_loops()
        self.loop_class = self.__classify_by_loops()

    def __classify_onzm(self) -> Optional[ONZM]:
        if len(self.tetrads) == 1:
            return None
        if any([t.onz is None for t in self.tetrads]):
            return None

        counter = Counter([t.onz.value[0] for t in self.tetrads])
        onz, support = counter.most_common()[0]
        if support != len(self.tetrads):
            onz = 'M'

        counter = Counter([tp.direction.value[0] for tp in self.tetrad_pairs])
        direction, support = counter.most_common()[0]
        if support != len(self.tetrad_pairs):
            direction = 'h'

        counter = Counter([t.onz.value[1] for t in self.tetrads])
        plus_minus, support = counter.most_common()[0]
        if support != len(self.tetrads):
            plus_minus = '*'

        return ONZM.from_value(f'{onz}{direction}{plus_minus}')

    def __classify_by_gba(self) -> List[GbaQuadruplexClassification]:
        gbas = set()
        for t in self.tetrads:
            gba = t.gba_class
            if gba is not None:
                gbas.add(gba.value[:-1])  # discard 'a' or 'b' subvariant
        roman_numerals = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8}
        gbas = sorted(gbas, key=lambda gba: roman_numerals.get(gba, 100))
        return list(map(lambda x: GbaQuadruplexClassification[x], gbas))

    def __find_tracts(self) -> List[Tract]:
        tracts = [[self.tetrads[0].nt1], [self.tetrads[0].nt2], [self.tetrads[0].nt3], [self.tetrads[0].nt4]]
        if len(self.tetrad_pairs) > 0:
            for tetrad_pair in self.tetrad_pairs:
                nt_dict = {
                    tetrad_pair.tetrad1.nt1: tetrad_pair.tetrad2_nts_best_order[0],
                    tetrad_pair.tetrad1.nt2: tetrad_pair.tetrad2_nts_best_order[1],
                    tetrad_pair.tetrad1.nt3: tetrad_pair.tetrad2_nts_best_order[2],
                    tetrad_pair.tetrad1.nt4: tetrad_pair.tetrad2_nts_best_order[3],
                }
                for i in range(4):
                    tracts[i].append(nt_dict[tracts[i][-1]])
        return [Tract(nts) for nts in tracts]

    def __find_loops(self) -> List[Loop]:
        if len(self.tetrads) == 1:
            return []

        loops = []
        tetrad_nucleotides = sorted([nt for tetrad in self.tetrads for nt in tetrad.nucleotides],
                                    key=lambda nt: nt.index)

        for i in range(1, len(tetrad_nucleotides)):
            nprev = tetrad_nucleotides[i - 1]
            ncur = tetrad_nucleotides[i]
            if ncur.index - nprev.index > 1 and ncur.chain == nprev.chain:
                for tract in self.tracts:
                    if nprev in tract.nucleotides and ncur in tract.nucleotides:
                        break
                else:
                    nts = list(filter(lambda nt: nprev.index < nt.index < ncur.index, self.structure3d.residues))
                    loop_type = self.__detect_loop_type(nprev, ncur)
                    loops.append(Loop(nts, loop_type))
        return loops

    def __detect_loop_type(self, nt_first: Residue3D, nt_last: Residue3D) -> Optional[LoopType]:
        tetrad_with_first = self.__find_tetrad_with_nt(nt_first)
        tetrad_with_last = self.__find_tetrad_with_nt(nt_last)

        if tetrad_with_first is None or tetrad_with_last is None:
            logging.warning(f'Failed to classify the loop between {nt_first} and {nt_last}')
            return None

        if tetrad_with_first == tetrad_with_last:
            # diagonal or laterals happen when first and last nt of a loop is in the same tetrad
            sign = self.__detect_loop_sign(nt_first, nt_last, tetrad_with_first)
            if sign is not None:
                return LoopType.from_value(f'lateral{sign}')
            return LoopType.diagonal

        tract_with_last = self.__find_tract_with_nt(nt_last)
        if tract_with_last is not None:
            # search along the tract to check what pairs with nt_first
            for nt in tract_with_last.nucleotides:
                if nt in tetrad_with_first.nucleotides:
                    sign = self.__detect_loop_sign(nt_first, nt, tetrad_with_first)
                    if sign is not None:
                        return LoopType.from_value(f'propeller{sign}')
        logging.warning(f'Failed to classify the loop between {nt_first} and {nt_last}')
        return None

    def __find_tetrad_with_nt(self, nt: Residue3D) -> Optional[Tetrad]:
        for tetrad in self.tetrads:
            if nt in tetrad.nucleotides:
                return tetrad
        return None

    def __find_tract_with_nt(self, nt: Residue3D) -> Optional[Tract]:
        for tract in self.tracts:
            if nt in tract.nucleotides:
                return tract
        return None

    def __detect_loop_sign(self, first: Residue3D, last: Residue3D, tetrad: Tetrad) -> Optional[str]:
        for pair in [tetrad.pair_12, tetrad.pair_23, tetrad.pair_34, tetrad.pair_41]:
            # main check
            if pair.nt1 == first and pair.nt2 == last:
                if pair.score() < pair.reverse().score():
                    return '-'
                return '+'
            # reverse check
            if pair.nt1 == last and pair.nt2 == first:
                if pair.score() < pair.reverse().score():
                    return '+'
                return '-'
        return None

    def __classify_by_loops(self) -> Optional[LoopClassification]:
        if len(self.loops) != 3 or any([loop.loop_type is None for loop in self.loops]):
            return None

        loop_classes = {
            'ppp': '1',
            'ppl': '2',
            'plp': '3',
            'lpp': '4',
            'pdp': '5',
            'lll': '6',
            'llp': '7',
            'lpl': '8',
            'pll': '9',
            'pdl': '10',
            'ldl': '11',
            'dpd': '12',
            'ldp': '13'
        }
        fingerprint = ''.join([loop.loop_type.value[0] for loop in self.loops])
        if fingerprint not in loop_classes:
            logging.error(f'Unknown loop classification: {fingerprint}')
            return None
        subtype = 'a' if self.loops[0 if fingerprint != 'dpd' else 1].loop_type.value[-1] == '-' else 'b'
        return LoopClassification.from_value(f'{loop_classes[fingerprint]}{subtype}')

    def __str__(self):
        builder = ''
        if len(self.tetrads) == 1:
            builder += '  single tetrad\n'
            builder += str(self.tetrads[0])
        else:
            builder += f'  {self.onzm.value if self.onzm is not None else "R"}'
            builder += f' {",".join(map(lambda gba: gba.value, self.gba_classes))}'
            if self.loop_class:
                builder += f' {self.loop_class.value} {self.loop_class.loop_progression()}'
            else:
                builder += f' n/a'
            builder += f' quadruplex with {len(self.tetrads)} tetrads\n'
            builder += str(self.tetrad_pairs[0].tetrad1)
            for tetrad_pair in self.tetrad_pairs:
                builder += str(tetrad_pair)
                builder += str(tetrad_pair.tetrad2)
            if self.tracts:
                builder += '\n    Tracts:\n'
                for tract in self.tracts:
                    builder += f'{tract}\n'
            if self.loops:
                builder += '\n    Loops:\n'
                for loop in self.loops:
                    builder += f'{loop}\n'
            builder += '\n'
        return builder


@dataclass
class Helix:
    tetrads: List[Tetrad]
    tetrad_pairs: List[TetradPair]
    structure3d: Structure3D
    quadruplexes: List[Quadruplex] = field(init=False)

    def __post_init__(self):
        self.quadruplexes = self.__find_quadruplexes()

    def __find_quadruplexes(self):
        if len(self.tetrad_pairs) == 0:
            return [Quadruplex(self.tetrads, [], self.structure3d)]

        quadruplexes = list()
        tetrads = list()

        for tetrad in [self.tetrad_pairs[0].tetrad1] + [tetrad_pair.tetrad2 for tetrad_pair in self.tetrad_pairs]:
            if tetrads:
                if tetrad.chains().isdisjoint(tetrads[-1].chains()):
                    quadruplexes.append(Quadruplex(tetrads, self.__filter_tetrad_pairs(tetrads), self.structure3d))
                    tetrads = list()
            tetrads.append(tetrad)

        quadruplexes.append(Quadruplex(tetrads, self.__filter_tetrad_pairs(tetrads), self.structure3d))

        return quadruplexes

    def __filter_tetrad_pairs(self, tetrads: List[Tetrad]) -> List[TetradPair]:
        chains = set()
        for tetrad in tetrads:
            chains.update(tetrad.chains())

        def check_tetrad(t: Tetrad) -> bool:
            return not t.chains().isdisjoint(chains)

        def check_pair(tp: TetradPair) -> bool:
            return check_tetrad(tp.tetrad1) and check_tetrad(tp.tetrad2)

        return list(filter(check_pair, self.tetrad_pairs))

    def __str__(self):
        builder = ''
        if len(self.tetrads) > 1:
            builder += f'n4-helix with {len(self.tetrads)} tetrads\n'
            for quadruplex in self.quadruplexes:
                builder += str(quadruplex)
        elif len(self.tetrads) == 1:
            builder += 'single tetrad without stacking\n'
            builder += str(self.tetrads[0])
        return builder


@dataclass
class Analysis:
    structure2d: Structure2D
    structure3d: Structure3D
    strict: bool
    no_reorder: bool
    stacking_mismatch: int
    base_pairs: List[BasePair3D] = field(init=False)
    base_pair_graph: Dict[Residue3D, List[Residue3D]] = field(init=False)
    base_pair_dict: Dict[Tuple[Residue3D, Residue3D], BasePair3D] = field(init=False)
    stacking_graph: Dict[Residue3D, List[Residue3D]] = field(init=False)
    tetrads: List[Tetrad] = field(init=False)
    tetrad_scores: Dict[Tetrad, Dict[Tetrad, Tuple[int, Tuple, Tuple]]] = field(init=False)
    tetrad_pairs: List[TetradPair] = field(init=False)
    helices: List[Helix] = field(init=False)
    ions: List[Atom3D] = field(init=False)
    sequence: str = field(init=False)
    line1: str = field(init=False)
    line2: str = field(init=False)
    shifts: Dict[Residue3D, int] = field(init=False)

    def __post_init__(self):
        self.base_pairs = self.structure3d.base_pairs(self.structure2d)
        self.base_pair_graph = self.structure3d.base_pair_graph(self.structure2d, self.strict)
        self.base_pair_dict = self.structure3d.base_pair_dict(self.structure2d, self.strict)
        self.stacking_graph = self.structure3d.stacking_graph(self.structure2d)
        self.tetrads = self.__find_tetrads(self.no_reorder)
        self.tetrad_scores = self.__calculate_tetrad_scores()
        self.tetrad_pairs = self.__find_tetrad_pairs(self.stacking_mismatch)
        self.helices = self.__find_helices()

        if not self.no_reorder:
            self.__find_best_chain_order()

        self.sequence, self.line1, self.line2, self.shifts = self.__generate_twoline_dotbracket()
        self.ions = self.__find_ions()
        self.__assign_ions_to_tetrads()

    def __find_tetrads(self, no_reorder=False) -> List[Tetrad]:
        # search for a tetrad: i -> j -> k -> l
        #                      ^--------------^
        tetrads = []
        for i in self.base_pair_graph:
            for j in filter(lambda x: x != i, self.base_pair_graph[i]):
                for k in filter(lambda x: x not in (i, j), self.base_pair_graph[j]):
                    for l in filter(lambda x: x not in (i, j, k) and i in self.base_pair_graph[x],
                                    self.base_pair_graph[k]):
                        if Tetrad.is_valid(i, j, k, l, self.base_pair_dict):
                            pair_12 = self.base_pair_dict[(i, j)]
                            pair_23 = self.base_pair_dict[(j, k)]
                            pair_34 = self.base_pair_dict[(k, l)]
                            pair_41 = self.base_pair_dict[(l, i)]
                            tetrads.append(Tetrad(i, j, k, l, pair_12, pair_23, pair_34, pair_41))

        # build graph of tetrads
        while tetrads:
            graph = defaultdict(list)
            for (ti, tj) in itertools.combinations(tetrads, 2):
                if not ti.is_disjoint(tj):
                    graph[ti].append(tj)
                    graph[tj].append(ti)

            # remove tetrad which conflicts the most with others
            # in case of a tie, remove one which has the worst planarity deviation
            candidates = sorted(tetrads, key=lambda t: (len(graph[t]), t.planarity_deviation),
                                reverse=True)
            if len(graph[candidates[0]]) > 0:
                tetrads.remove(candidates[0])
            else:
                break

        return sorted(tetrads, key=lambda t: min(map(lambda nt: nt.index, t.nucleotides)))

    def __calculate_tetrad_scores(self) \
            -> Dict[Tetrad, Dict[Tetrad, Tuple[int, Tuple, Tuple]]]:
        def is_next_by_stacking(nt1: Residue3D, nt2: Residue3D) -> bool:
            return nt2 in self.stacking_graph.get(nt1, [])

        def is_next_sequentially(nt1: Residue3D, nt2: Residue3D) -> bool:
            return nt1.chain == nt2.chain and abs(nt1.index - nt2.index) == 1

        tetrad_scores = defaultdict(dict)

        for ti, tj in itertools.combinations(self.tetrads, 2):
            nts1 = ti.nucleotides
            best_score = 0
            best_score_sequential = 0
            best_score_stacking = 0
            best_order = tj.nucleotides

            n1, n2, n3, n4 = tj.nucleotides
            viable_permutations = [(n1, n2, n3, n4), (n2, n3, n4, n1), (n3, n4, n1, n2), (n4, n1, n2, n3),
                                   (n1, n4, n3, n2), (n4, n3, n2, n1), (n3, n2, n1, n4), (n2, n1, n4, n3)]

            for nts2 in viable_permutations:
                score_stacking = [1 if is_next_by_stacking(nts1[i], nts2[i]) else 0 for i in range(4)]
                score_sequential = [1 if is_next_sequentially(nts1[i], nts2[i]) else 0 for i in range(4)]
                score = sum([max(score_stacking[i], score_sequential[i]) for i in range(4)])
                score_sequential = sum(score_sequential)
                score_stacking = sum(score_stacking)

                if (score, score_sequential, score_stacking) > (best_score, best_score_sequential, best_score_stacking):
                    best_score, best_score_sequential, best_score_stacking = score, score_sequential, score_stacking
                    best_order = nts2
                if best_score == 4:
                    break

            tetrad_scores[ti][tj] = (best_score, nts1, best_order)
            tetrad_scores[tj][ti] = (best_score, best_order, nts1)

        return tetrad_scores

    def __find_tetrad_pairs(self, stacking_mismatch: int) -> List[TetradPair]:
        tetrads = list(self.tetrads)
        best_score = 0
        best_order = tetrads

        for ti in tetrads:
            score = 0
            order = [ti]
            candidates = set(self.tetrads) - {ti}

            while candidates:
                tj = max([tj for tj in candidates], key=lambda tk: self.tetrad_scores[ti][tk][0])
                score += self.tetrad_scores[ti][tj][0]
                order.append(tj)
                candidates.remove(tj)
                ti = tj

            if score > best_score:
                best_score = score
                best_order = order

            if best_score == (len(self.tetrads) - 1) * 4:
                break

        tetrad_pairs = []

        for i in range(1, len(best_order)):
            ti, tj = best_order[i - 1], best_order[i]
            score = self.tetrad_scores[ti][tj][0]

            if score >= (4 - stacking_mismatch):
                nts1, nts2 = self.tetrad_scores[ti][tj][1:]
                stacked = {nts1[i]: nts2[i] for i in range(4)}
                stacked.update({v: k for k, v in stacked.items()})
                tetrad_pairs.append(TetradPair(ti, tj, stacked))
                order = (stacked[ti.nt1], stacked[ti.nt2], stacked[ti.nt3], stacked[ti.nt4])
                tj.reorder_to_match_other_tetrad(order)

        return tetrad_pairs

    def __find_helices(self):
        helices = []
        helix_tetrads = []
        helix_tetrad_pairs = []

        for tp in self.tetrad_pairs:
            ti, tj = tp.tetrad1, tp.tetrad2
            if not helix_tetrads:
                helix_tetrads.append(ti)
            score = self.tetrad_scores[helix_tetrads[-1]][tj][0]
            if score >= (4 - self.stacking_mismatch):
                helix_tetrads.append(tj)
                helix_tetrad_pairs.append(tp)
            else:
                helices.append(Helix(helix_tetrads, helix_tetrad_pairs, self.structure3d))
                helix_tetrads = [ti, tj]
                helix_tetrad_pairs = [tp]

        if helix_tetrads:
            helices.append(Helix(helix_tetrads, helix_tetrad_pairs, self.structure3d))

        for tetrad in self.tetrads:
            if not any([tetrad in helix.tetrads for helix in helices]):
                helices.append(Helix([tetrad], [], self.structure3d))

        return helices

    def __find_best_chain_order(self):
        chain_groups = self.__group_related_chains()
        final_order = []

        for chains in chain_groups:
            best_permutation, best_score = chains, (1e10, 1e10)

            if len(chains) > 1:
                for permutation in itertools.permutations(chains):
                    self.__reorder_chains(permutation)
                    classifications = [t.onz for h in self.helices for t in h.tetrads]
                    logging.debug(
                        f'Checking reorder: {" ".join(permutation)} {" ".join(map(lambda c: c.value, classifications))}')

                    onz_score = sum(c.score() for c in classifications)
                    chain_order_score = self.__chain_order_score(permutation)
                    score = (onz_score, chain_order_score)

                    if score < best_score:
                        best_score = score
                        best_permutation = permutation
                    elif score == best_score:
                        # in case of a tie, pick permutation earlier in lexicographical sense
                        if permutation < best_permutation:
                            best_permutation = permutation

            final_order.extend(best_permutation)

        if len(final_order) > 1:
            self.__reorder_chains(final_order)
            classifications = [t.onz for h in self.helices for t in h.tetrads]
            logging.debug(f'Selected chain order: {" ".join(final_order)} '
                          f'{" ".join(map(lambda onz: onz.value, classifications))}')

            self.tetrads = self.__find_tetrads(True)
            self.tetrad_scores = self.__calculate_tetrad_scores()
            self.tetrad_pairs = self.__find_tetrad_pairs(self.stacking_mismatch)
            self.helices = self.__find_helices()

    def __group_related_chains(self) -> List[List[str]]:
        candidates = set()

        for h in self.helices:
            for t in h.tetrads:
                candidates.add(frozenset([t.nt1.chain, t.nt2.chain, t.nt3.chain, t.nt4.chain]))

        candidates = [set(c) for c in candidates]
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

        candidates = sorted(candidates, key=lambda x: len(x), reverse=True)
        groups = []

        for candidate in candidates:
            if any([group.issuperset(candidate) for group in groups]):
                continue
            groups.append(candidate)

        return sorted([sorted(group) for group in groups], key=lambda x: x[0])

    def __reorder_chains(self, chain_order: Iterable[str]):
        i = 1
        for chain in chain_order:
            for nt in self.structure3d.residues:
                if nt.chain == chain:
                    nt.index = i
                    i += 1
        for nt in self.structure3d.residues:
            if nt.chain not in chain_order:
                nt.index = i
                i += 1

        if len(self.tetrad_pairs) > 0:
            self.tetrad_pairs[0].tetrad1.reorder_to_match_5p_3p()
            for tp in self.tetrad_pairs:
                order = (tp.stacked[tp.tetrad1.nt1], tp.stacked[tp.tetrad1.nt2],
                         tp.stacked[tp.tetrad1.nt3], tp.stacked[tp.tetrad1.nt4])
                tp.tetrad2.reorder_to_match_5p_3p()  # this is required to recalculate ONZ
                tp.tetrad2.reorder_to_match_other_tetrad(order)

    def __chain_order_score(self, chain_order: Tuple[str, ...]) -> int:
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

    def __find_ions(self) -> List[Atom3D]:
        metal_atom_names = set([ion.value.upper() for ion in Ion])
        ions = []
        used = set()
        for residue in self.structure3d.residues:
            for atom in residue.atoms:
                if atom.atomName.upper() in metal_atom_names:
                    coordinates = tuple(atom.coordinates())
                    if coordinates not in used:
                        ions.append(atom)
                        used.add(coordinates)
        return ions

    def __assign_ions_to_tetrads(self) \
            -> Tuple[Dict[Tetrad, List[Atom3D]], Dict[Tuple[Tetrad, Residue3D], List[Atom3D]]]:
        if len(self.tetrads) == 0:
            return {}, {}

        ions_channel = defaultdict(list)
        ions_outside = defaultdict(list)

        for ion in self.ions:
            min_distance = math.inf
            min_tetrad = self.tetrads[0]

            for tetrad in self.tetrads:
                distance = numpy.linalg.norm(ion.coordinates() - tetrad.center())
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
                        distance = numpy.linalg.norm(ion.coordinates() - atom.coordinates())
                        if distance < min_distance:
                            min_distance = distance
                            min_tetrad = tetrad
                            min_nt = nt

            # TODO: verify threshold of 3A between an ion and an atom
            if min_distance < 3.0:
                ions_outside[(min_tetrad, min_nt)].append(ion)
                continue

            logging.debug(f'Skipping an ion, because it is too far from any tetrad (distance={min_distance})')

        for tetrad, ions in ions_channel.items():
            tetrad.ions_channel = ions
        for pair, ions in ions_outside.items():
            tetrad, residue = pair
            tetrad.ions_outside[residue] = ions

    def __generate_twoline_dotbracket(self) -> Tuple[str, str, str, Dict[Residue3D, int]]:
        layer1, layer2 = [], []
        for tetrad in self.tetrads:
            layer1.extend([tetrad.pair_12, tetrad.pair_34])
            layer2.extend([tetrad.pair_23, tetrad.pair_41])
        sequence, line1, shifts = self.__elimination_conflicts(layer1)
        _, line2, _ = self.__elimination_conflicts(layer2)
        return sequence, line1, line2, shifts

    def __elimination_conflicts(self, pairs: List[BasePair3D]) -> Tuple[str, str, Dict[Residue3D, int]]:
        orders = dict()
        order = 0
        queue = list(pairs)
        removed = []

        while queue:
            conflicts = defaultdict(list)
            for pi, pj in itertools.combinations(queue, 2):
                if pi.conflicts_with(pj):
                    conflicts[pi].append(pj)
                    conflicts[pj].append(pi)
            if conflicts:
                pair, _ = max(conflicts.items(), key=lambda x: (len(x[1]), x[0].nt1))
                removed.append(pair)
                queue.remove(pair)
            else:
                orders.update({pair: order for pair in queue})
                queue, removed = removed, []
                order += 1

        opening = '([{<' + string.ascii_uppercase
        closing = ')]}>' + string.ascii_lowercase
        dotbracket = dict()
        for pair, order in orders.items():
            nt1, nt2 = sorted([pair.nt1, pair.nt2])
            dotbracket[nt1] = opening[order]
            dotbracket[nt2] = closing[order]

        sequence = ''
        structure = ''
        shifts = dict()
        shift_value = 0
        chain = None
        for nt in sorted(filter(lambda nt: nt.is_nucleotide, self.structure3d.residues), key=lambda nt: nt.index):
            if chain and chain != nt.chain:
                sequence += '-'
                structure += '-'
                shift_value += 1
            sequence += nt.one_letter_name
            structure += dotbracket.get(nt, '.')
            shifts[nt] = shift_value
            chain = nt.chain
        return sequence, structure, shifts

    def __str__(self):
        builder = f'Chain order: {" ".join(self.__chain_order())}\n'
        for helix in self.helices:
            builder += str(helix)
        builder += f'{self.sequence}\n{self.line1}\n{self.line2}'
        return builder

    def __chain_order(self) -> List[str]:
        only_nucleic_acids = filter(lambda nt: nt.is_nucleotide, self.structure3d.residues)
        return list({nt.chain: 0 for nt in sorted(only_nucleic_acids, key=lambda nt: nt.index)}.keys())

    def canonical(self) -> List[BasePair3D]:
        return [base_pair for base_pair in self.base_pairs if base_pair.is_canonical()]


@dataclass
class Visualizer:
    analysis: Analysis
    tetrads: List[Tetrad]
    complete2d: bool
    onz_dict: Dict[BasePair3D, ONZ] = field(init=False)

    def __post_init__(self):
        self.onz_dict = {pair: tetrad.onz for tetrad in self.tetrads for pair in
                         [tetrad.pair_12, tetrad.pair_23, tetrad.pair_34, tetrad.pair_41]}

    def visualize(self, prefix: str, suffix: str):
        fasta = tempfile.NamedTemporaryFile('w+', suffix='.fasta')
        fasta.write(f'>{prefix}-{suffix}\n')
        fasta.write(self.analysis.sequence)
        fasta.flush()

        layer1, layer2 = [], []
        for tetrad in self.tetrads:
            layer1.extend([tetrad.pair_12, tetrad.pair_34])
            layer2.extend([tetrad.pair_23, tetrad.pair_41])
        helix1 = self.__to_helix(layer1, self.analysis.canonical() if self.complete2d else [])
        helix2 = self.__to_helix(layer2)

        currdir = os.path.dirname(os.path.realpath(__file__))
        output_pdf = f'{prefix}-{suffix}.pdf'
        run = subprocess.run([os.path.join(currdir, 'quadraw.R'), fasta.name, helix1.name, helix2.name, output_pdf],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if run.returncode == 0:
            print('\nPlot:', output_pdf)
        else:
            logging.error(f'Failed to prepare visualization, reason:\n  {run.stderr.decode()}')

    def __to_helix(self, layer: List[BasePair3D],
                   canonical: Optional[List[BasePair3D]] = None) -> tempfile.NamedTemporaryFile():
        onz_value = {ONZ.O_PLUS: 1, ONZ.O_MINUS: 2, ONZ.N_PLUS: 3, ONZ.N_MINUS: 4, ONZ.Z_PLUS: 5, ONZ.Z_MINUS: 6}
        nucleotides = self.analysis.structure3d.residues
        shifts = self.analysis.shifts

        helix = tempfile.NamedTemporaryFile('w+', suffix='.helix')
        helix.write(f'#{len(self.analysis.sequence) + 1}\n')
        helix.write('i\tj\tlength\tvalue\n')

        for pair in layer:
            x, y = pair.nt1, pair.nt2
            x, y = nucleotides.index(x) + 1 + shifts[x], nucleotides.index(y) + 1 + shifts[y]
            onz = self.onz_dict[pair]
            helix.write(f'{x}\t{y}\t1\t{onz_value.get(onz, 7)}\n')
        if canonical:
            for pair in canonical:
                x, y = pair.nt1, pair.nt2
                x, y = nucleotides.index(x) + 1 + shifts[x], nucleotides.index(y) + 1 + shifts[y]
                helix.write(f'{x}\t{y}\t1\t8\n')

        helix.flush()
        return helix


class AnalysisSimple:
    def __init__(self, structure2d: Structure2D, structure3d: Structure3D):
        self.pairs: List[BasePair3D] = structure3d.base_pairs(structure2d)
        self.graph: Dict[Residue3D, List[Residue3D]] = structure3d.base_pair_graph(structure2d)
        self.pair_dict: Dict[Tuple[Residue3D, Residue3D], BasePair3D] = structure3d.base_pair_dict(structure2d)

    def has_tetrads(self):
        tetrads = set()
        for i in self.graph:
            for j in filter(lambda x: x != i, self.graph[i]):
                for k in filter(lambda x: x not in (i, j), self.graph[j]):
                    for l in filter(lambda x: x not in (i, j, k) and x in self.graph[i], self.graph[k]):
                        if Tetrad.is_valid(i, j, k, l, self.pair_dict):
                            tetrads.add(frozenset([i, j, k, l]))
                        if len(tetrads) > 1:
                            return True
        return False


def center_of_mass(atoms):
    coords = [atom.coordinates() for atom in atoms]
    xs = (coord[0] for coord in coords)
    ys = (coord[1] for coord in coords)
    zs = (coord[2] for coord in coords)
    return numpy.array((sum(xs) / len(coords), sum(ys) / len(coords), sum(zs) / len(coords)))


def eltetrado(structure2d: Structure2D, structure3d: Structure3D, strict: bool, no_reorder: bool,
              stacking_mismatch: int) -> Analysis:
    return Analysis(structure2d, structure3d, strict, no_reorder, stacking_mismatch)


def has_tetrad(structure2d: Structure2D, structure3d: Structure3D) -> bool:
    structure = AnalysisSimple(structure2d, structure3d)
    return structure.has_tetrads()
