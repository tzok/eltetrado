import itertools
import logging
import math
import numpy
import os
import string
import subprocess
import tempfile

from collections import defaultdict, Counter
from typing import Dict, Iterable, List, Tuple, FrozenSet, Set, Optional
from eltetrado.compute.model import Atom3D, Structure3D, Residue3D

logging.basicConfig(level=os.environ.get("LOGLEVEL", "INFO"))

LW_SCORE = {
    'cWW': 1,
    'tWW': 2,
    'cWH': 3,
    'tWH': 4,
    'cWS': 5,
    'tWS': 6,
    'cHW': 7,
    'tHW': 8,
    'cHH': 9,
    'tHH': 10,
    'cHS': 11,
    'tHS': 12,
    'cSW': 13,
    'tSW': 14,
    'cSH': 15,
    'tSH': 16,
    'cSS': 17,
    'tSS': 18
}


class Nucleotide:
    """
    Metadata about a single nucleotide.

    Attributes:
        outermost_atoms     Dict representing expected name of atom involved in glycosidic bond
        innermost_atoms     Dist representing expected name of atom closest to the tetrad center
    """
    outermost_atoms = {'A': 'N9', 'G': 'N9', 'C': 'N1', 'U': 'N1', 'T': 'N1'}
    innermost_atoms = {'A': 'N6', 'G': 'O6', 'C': 'N4', 'U': 'O4', 'T': 'O4'}

    @staticmethod
    def detect_molecule(resname: str):
        if resname in ('DA', 'DC', 'DG', 'DT'):
            return 'DNA'
        if resname in ('A', 'C', 'G', 'U'):
            return 'RNA'
        return 'Other'

    @staticmethod
    def syn_anti(chi: float) -> str:
        if chi is None:
            return '...'
        return 'syn' if -90 < chi < 90 else 'anti'

    def __init__(self, nt: dict, structure3d: Structure3D):
        self.model: int = 1 if nt['nt_id'].find(':') == -1 else int(nt['nt_id'].split(':')[0])
        self.chain: str = nt['chain_name']
        self.number: int = nt['nt_resnum']
        self.icode: str = ' ' if nt['nt_id'].find('^') == -1 else nt['nt_id'].split('^')[1]
        self.molecule: str = self.detect_molecule(nt['nt_name'])
        self.full_name: str = nt['nt_id']
        self.short_name: str = nt['nt_code']
        self.chi: float = nt['chi']
        self.glycosidic_bond: str = self.syn_anti(self.chi)
        self.index: int = nt['index']
        self.model_chain: str = self.chain if nt['nt_id'].find(':') == -1 else f'{self.model}:{self.chain}'
        self.residue3d: Optional[Residue3D] = None

        if structure3d:
            for residue in structure3d.residues:
                if residue.chain_identifier == self.chain \
                        and residue.residue_number == self.number \
                        and residue.insertion_code == self.icode:
                    self.residue3d = residue
                    break
            else:
                logging.error(f'Failed to find data in PDB or mmCIF file for residue {self.full_name}')

    def __eq__(self, other):
        return self.full_name == other.full_name

    def __hash__(self):
        return hash(self.full_name)

    def __lt__(self, other):
        return self.index < other.index

    def __repr__(self):
        return self.full_name

    def __str__(self):
        return self.full_name

    def find_atom(self, expected_atom: str) -> Optional[Atom3D]:
        if self.residue3d:
            for atom in self.residue3d.atoms:
                if atom.atom_name == expected_atom:
                    return atom
        return None

    def is_purine(self):
        return self.find_atom('N9') is not None

    def __outer_generator(self):
        # try to find expected atom name
        upper = self.short_name.upper()
        if upper in self.outermost_atoms:
            yield self.find_atom(self.outermost_atoms[upper])

        # try to get generic name for purine/pyrimidine
        if self.is_purine():
            yield self.find_atom('N9')
        else:
            yield self.find_atom('N1')

        # try to find at least C1' next to nucleobase
        yield self.find_atom("C1'")

        # get any atom
        if self.residue3d.atoms:
            yield self.residue3d.atoms[0]

        # last resort, create pseudoatom at (0, 0, 0)
        logging.error(
            f'Failed to determine the outermost atom for nucleotide {self}, so an arbitrary atom will be used')
        yield Atom3D('', '', '', 1, ' ', 0.0, 0.0, 0.0)

    def __inner_generator(self):
        # try to find expected atom name
        upper = self.short_name.upper()
        if upper in self.innermost_atoms:
            yield self.find_atom(self.innermost_atoms[upper])

        # try to get generic name for purine/pyrimidine
        if self.is_purine():
            yield self.find_atom('C6')
        else:
            yield self.find_atom('C4')

        # try to find any atom at position 4 or 6 for purine/pyrimidine respectively
        if self.is_purine():
            yield self.find_atom('O6')
            yield self.find_atom('N6')
            yield self.find_atom('S6')
        else:
            yield self.find_atom('O4')
            yield self.find_atom('N4')
            yield self.find_atom('S4')

        # get any atom
        if self.residue3d.atoms:
            yield self.residue3d.atoms[0]

        # last resort, create pseudoatom at (0, 0, 0)
        logging.error(
            f'Failed to determine the innermost atom for nucleotide {self}, so an arbitrary atom will be used')
        yield Atom3D('', '', '', 1, ' ', 0.0, 0.0, 0.0)

    def outermost_atom(self) -> Atom3D:
        return next(filter(None, self.__outer_generator()))

    def innermost_atom(self) -> Atom3D:
        return next(filter(None, self.__inner_generator()))


class Pair:

    @staticmethod
    def is_valid(lw: str):
        return len(lw) == 3 and lw[0] in 'ct' and lw[1] in 'WHS' and lw[2] in 'WHS'

    def __init__(self, nt1: Nucleotide, nt2: Nucleotide, lw: str, saenger: str):
        self.pair: Tuple[Nucleotide, Nucleotide] = (nt1, nt2)
        self.lw: str = lw
        self.saenger: str = saenger

    def __eq__(self, other):
        return self.pair == other.pair and self.lw == other.lw and self.saenger == other.saenger

    def __hash__(self):
        return hash((self.pair, self.lw, self.saenger))

    def __str__(self):
        return f'{self.pair[0]} {self.pair[1]} {self.lw} {self.score()}'

    def conflicts_with(self, other) -> bool:
        xi, yi = sorted((self.pair[0].index, self.pair[1].index))
        xj, yj = sorted((other.pair[0].index, other.pair[1].index))
        return xi < xj < yi < yj or xj < xi < yj < yi

    def reverse(self):
        lw = f'{self.lw[0]}{self.lw[2]}{self.lw[1]}'
        return Pair(self.pair[1], self.pair[0], lw, self.saenger)

    def score(self) -> int:
        return LW_SCORE.get(self.lw, 100)

    def is_canonical(self) -> bool:
        return self.saenger in ('19-XIX', '20-XX', '28-XXVIII')


class Tetrad:

    @staticmethod
    def is_valid(nt1: Nucleotide, nt2: Nucleotide, nt3: Nucleotide, nt4: Nucleotide,
                 pairs: Dict[Tuple[Nucleotide, Nucleotide], Pair]) -> bool:
        p1, p2, p3, p4 = pairs[(nt1, nt2)], pairs[(nt2, nt3)], pairs[(nt3, nt4)], pairs[(nt4, nt1)]
        for pi, pj in ((p1, p4), (p2, p1), (p3, p2), (p4, p3)):
            if pi.lw[1] == pj.lw[2]:
                return False
        return True

    def __init__(self,
                 nt1: Nucleotide,
                 nt2: Nucleotide,
                 nt3: Nucleotide,
                 nt4: Nucleotide,
                 pairs: Dict[Tuple[Nucleotide, Nucleotide], Pair],
                 no_reorder=False):
        self.nucleotides: Tuple[Nucleotide, Nucleotide, Nucleotide, Nucleotide] = (nt1, nt2, nt3, nt4)
        self.pairs: Tuple[Pair, Pair, Pair,
                          Pair] = (pairs[(nt1, nt2)], pairs[(nt2, nt3)], pairs[(nt3, nt4)], pairs[(nt4, nt1)])
        self.set: FrozenSet[Nucleotide] = frozenset(self.nucleotides)
        self.planarity_deviation: float = self.__calculate_planarity_deviation()
        self.no_reorder = no_reorder
        self.chains: Set[str] = self.__chains()
        self.center: numpy.ndarray = self.__calculate_center_point()
        self.ions_channel: List[Atom3D] = []
        self.ions_outside: Dict[Nucleotide, List] = {}
        self.__score = sum(x.score() for x in self.pairs)
        self.__hash = hash(self.set)

    def __hash__(self):
        return self.__hash

    def __eq__(self, other):
        return hash(self) == hash(other)

    def __iter__(self):
        return iter(self.nucleotides)

    def __repr__(self):
        return '{}-{}-{}-{}'.format(*self.nucleotides)

    def __str__(self):
        return '    {} {} {} {} {}-{}-{}-{} {} {} ' \
               'planarity={} ' \
               'ions_channel={} ' \
               'ions_outside={}\n'.format(self.nucleotides[0], self.nucleotides[1], self.nucleotides[2],
                                          self.nucleotides[3], self.pairs[0].lw, self.pairs[1].lw, self.pairs[2].lw,
                                          self.pairs[3].lw, self.get_classification(), self.gba_classification(),
                                          round(self.planarity_deviation, 2),
                                          ','.join([atom.atom_name for atom in self.ions_channel]),
                                          {k: ','.join(map(lambda atom: atom.atom_name, v)) for k, v in
                                           self.ions_outside.items()})

    def stems_with(self, other) -> bool:
        index_diff = (j - i for i, j in zip(self.sorted_indices(), other.sorted_indices()))
        return all(abs(diff) == 1 for diff in index_diff)

    def count_non_stacked_bases(self, other, stacking: Set[Tuple[Nucleotide]]) -> Tuple[int, int]:
        ti, tj = set(self.nucleotides), set(other.nucleotides)

        for ni, nj in itertools.product(self.nucleotides, other.nucleotides):
            if ni not in ti or nj not in tj:
                continue
            for stack in stacking:
                if ni in stack and nj in stack and abs(stack.index(ni) - stack.index(nj)) == 1:
                    ti.remove(ni)
                    tj.remove(nj)
                    break

        return len(ti), len(tj)

    def sorted_indices(self) -> List[int]:
        return sorted(map(lambda x: x.index, self.nucleotides))

    def get_score(self) -> int:
        return self.__score

    def is_disjoint(self, other) -> bool:
        return self.set.isdisjoint(other.set)

    def reorder(self):
        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = (nt.index for nt in self.nucleotides)
        indices = sorted((ni, nj, nk, nl))
        ni, nj, nk, nl = (indices.index(x) for x in (ni, nj, nk, nl))

        nmin = min(ni, nj, nk, nl)
        if nmin == ni:
            pass
        elif nmin == nj:
            self.nucleotides = self.nucleotides[1:4] + self.nucleotides[0:1]
            self.pairs = self.pairs[1:4] + self.pairs[0:1]
        elif nmin == nk:
            self.nucleotides = self.nucleotides[2:4] + self.nucleotides[0:2]
            self.pairs = self.pairs[2:4] + self.pairs[0:2]
        else:
            self.nucleotides = self.nucleotides[3:4] + self.nucleotides[0:3]
            self.pairs = self.pairs[3:4] + self.pairs[0:3]

        # flip order if necessary
        if self.pairs[0].score() > self.pairs[3].reverse().score():
            self.nucleotides = self.nucleotides[0:1] + tuple(reversed(self.nucleotides[1:]))
            self.pairs = tuple(self.pairs[i].reverse() for i in (3, 2, 1, 0))

    def get_classification(self) -> str:
        if self.no_reorder and len(set((nt.chain for nt in self.nucleotides))) > 1:
            return 'n/a'

        self.reorder()

        # transform into (0, 1, 2, 3)
        ni, nj, nk, nl = (nt.index for nt in self.nucleotides)
        indices = sorted((ni, nj, nk, nl))
        ni, nj, nk, nl = (indices.index(x) for x in (ni, nj, nk, nl))

        order = (nj, nk, nl)
        if order == (1, 2, 3):
            return 'O+'
        elif order == (3, 2, 1):
            return 'O-'
        elif order == (1, 3, 2):
            return 'N+'
        elif order == (2, 3, 1):
            return 'N-'
        elif order == (2, 1, 3):
            return 'Z+'
        elif order == (3, 1, 2):
            return 'Z-'
        else:
            logging.error(f'Impossible combination: {ni} {nj} {nk} {nl}')
            return 'n/a'

    def gba_classification(self):
        """
        See: Webba da Silva, M. (2007). Geometric Formalism for DNA Quadruplex Folding.
        Chemistry - A European Journal, 13(35), 9738–9745. https://doi.org/10.1002/chem.200701255

        :return: Classification according to Webba da Silva or n/a
        """
        # without reordering, do not classify bi- and tetramolecular quadruplexes
        if self.no_reorder and len(set((nt.chain for nt in self.nucleotides))) > 1:
            return 'n/a'

        # without all nucleotides having a valid syn/anti, this classification is impossible
        if not all([nt.glycosidic_bond in ('syn', 'anti') for nt in self.nucleotides]):
            return 'n/a'

        # for Webba da Silva's classifcation, the original order needs to be maintained
        self.reorder()
        if self.nucleotides[1] < self.nucleotides[3]:
            nucleotides = self.nucleotides
        else:
            nucleotides = [self.nucleotides[0]] + list(reversed(self.nucleotides[1:]))

        # this will create a 4-letter string made of 's' for syn or 'a' for anti
        fingerprint = ''.join([nt.glycosidic_bond[0] for nt in nucleotides])

        # this dict has all classes mapped to fingerprints
        gba_classes = {
            'aass': 'Ia',
            'ssaa': 'Ib',
            'asas': 'IIa',
            'sasa': 'IIb',
            'asaa': 'IIIa',
            'sass': 'IIIb',
            'aaas': 'IVa',
            'sssa': 'IVb',
            'aasa': 'Va',
            'ssas': 'Vb',
            'assa': 'VIa',
            'saas': 'VIb',
            'asss': 'VIIa',
            'saaa': 'VIIb',
            'aaaa': 'VIIIa',
            'ssss': 'VIIIb'
        }

        if fingerprint not in gba_classes:
            logging.error(f'Impossible combination of syn/anti: {[nt.glycosidic_bond for nt in self.nucleotides]}')
            return 'n/a'
        return gba_classes[fingerprint]

    def __chains(self) -> Set[str]:
        return set([nt.chain for nt in self.nucleotides])

    def __calculate_planarity_deviation(self):
        if all(nt.residue3d for nt in self.nucleotides):
            outer = [nt.outermost_atom() for nt in self.nucleotides]
            inner = [nt.innermost_atom() for nt in self.nucleotides]
            if all(outer) and all(inner):
                return numpy.linalg.norm(center_of_mass(outer) - center_of_mass(inner))
        return float('nan')

    def __calculate_center_point(self):
        return sum(atom.coordinates() for atom in map(lambda nt: nt.innermost_atom(), self.nucleotides)) / len(
            self.nucleotides)


class TetradPair:

    def __init__(self, tetrad1: Tetrad, tetrad2: Tetrad, stacked: Dict[Nucleotide, Nucleotide]):
        self.tetrad1: Tetrad = tetrad1
        self.tetrad2: Tetrad = tetrad2
        self.stacked: Dict[Nucleotide, Nucleotide] = stacked
        self.direction: str = self.__determine_direction()
        self.rise: float = self.__calculate_rise()
        self.twist: float = self.__calculate_twist()

    def __str__(self):
        return f'      direction={self.direction} rise={round(self.rise, 2)} twist={round(self.twist, 2)}\n'

    def reverse(self):
        return TetradPair(self.tetrad2, self.tetrad1, self.stacked)

    def __determine_direction(self):
        # count directions 5' -> 3' as +1 or -1
        counter = Counter(
            1 if j - i > 0 else -1 for i, j in zip(self.tetrad1.sorted_indices(), self.tetrad2.sorted_indices()))
        direction, count = counter.most_common()[0]
        # all in the same direction
        if count == 4:
            return 'parallel'
        # two in +, one in - direction
        elif count == 2:
            return 'antiparallel'
        else:
            return 'hybrid'

    def __calculate_rise(self):
        f = Nucleotide.outermost_atom
        g = Nucleotide.innermost_atom
        mine = tuple(itertools.chain.from_iterable((f(nt), g(nt)) for nt in self.tetrad1.nucleotides))
        theirs = tuple(itertools.chain.from_iterable((f(nt), g(nt)) for nt in self.tetrad2.nucleotides))
        if all(mine) and all(theirs):
            return numpy.linalg.norm(center_of_mass(mine) - center_of_mass(theirs))
        return float('nan')

    def __calculate_twist(self):
        nt1_self = self.tetrad1.nucleotides[0]
        nt2_self = self.tetrad1.nucleotides[1]
        nt1_other = self.stacked[nt1_self]
        nt2_other = self.stacked[nt2_self]

        v1 = nt1_self.find_atom("C1'").coordinates() - nt2_self.find_atom("C1'").coordinates()
        v1 = v1 / numpy.linalg.norm(v1)
        v2 = nt1_other.find_atom("C1'").coordinates() - nt2_other.find_atom("C1'").coordinates()
        v2 = v2 / numpy.linalg.norm(v2)
        return math.degrees(numpy.arccos(numpy.clip(numpy.dot(v1, v2), -1.0, 1.0)))


class Loop:

    def __init__(self, nucleotides: List[Nucleotide], loop_type: str):
        self.nucleotides = nucleotides
        self.loop_type = loop_type

    def __iter__(self):
        return iter(self.nucleotides)

    def __str__(self):
        return f'      {self.loop_type} {", ".join(map(str, self.nucleotides))}'


class Tract:

    def __init__(self, nuleotides: List[Nucleotide]):
        self.nucleotides = nuleotides

    def __iter__(self):
        return iter(self.nucleotides)

    def __str__(self):
        return f'      {", ".join(map(lambda nt: nt.full_name, self.nucleotides))}'


class Quadruplex:

    def __init__(self, tetrads: List[Tetrad], tetrad_pairs: List[TetradPair], nucleotides: Dict[str, Nucleotide]):
        self.tetrads: List[Tetrad] = tetrads
        self.tetrad_pairs: List[TetradPair] = tetrad_pairs
        self.nucleotides: Dict[str, Nucleotide] = nucleotides

        self.tracts: List[Tract] = self.__tracts()
        self.loops: List[Loop] = self.__loops()

        self.gba_classification = self.__gba_classification()
        self.loop_classification = self.__loop_classification()

    def __str__(self):
        builder = ''
        if len(self.tetrads) == 1:
            builder += '  single tetrad\n'
            builder += str(self.tetrads[0])
        else:
            if any(t.get_classification() == 'n/a' for t in self.tetrads):
                builder += '  R {} {} quadruplex with {} tetrads\n'.format(','.join(self.gba_classification),
                                                                           self.loop_classification, len(self.tetrads))
            else:
                builder += '  {}{}{} {} {} quadruplex with {} tetrads\n'.format(self.onzm_classification(),
                                                                                self.direction(), self.sign(),
                                                                                ','.join(self.gba_classification),
                                                                                self.loop_classification,
                                                                                len(self.tetrads))
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

    def onzm_classification(self) -> str:
        classifications = [t.get_classification()[0] for t in self.tetrads]
        counter = Counter(classifications)
        onz, support = counter.most_common()[0]
        if support == len(self.tetrads):
            return onz[0]  # O, N or Z
        return 'M'

    def direction(self) -> str:
        if len(self.tetrads) == 1:
            return 'n/a'
        else:
            counter = Counter((pair.direction for pair in self.tetrad_pairs))
            direction, support = counter.most_common()[0]
            if support == len(self.tetrads) - 1:
                return direction[0]
            return 'h'

    def sign(self) -> str:
        signs = set((tetrad.get_classification()[1] for tetrad in self.tetrads))
        if len(signs - {'-', '+'}) > 0:
            logging.error(f'Tetrad classification different than [ONZ][+-]: '
                          f'{[tetrad.get_classification() for tetrad in self.tetrads]}')
            return '*'
        if len(signs) == 1:
            return signs.pop()
        return '*'

    def __gba_classification(self) -> List[str]:
        roman_numerals = {'I': 1, 'II': 2, 'III': 3, 'IV': 4, 'V': 5, 'VI': 6, 'VII': 7, 'VIII': 8}
        gbas = map(lambda tetrad: tetrad.gba_classification(), self.tetrads)
        gbas = filter(lambda gba: gba != 'n/a', gbas)
        gbas = map(lambda gba: gba[:-1], gbas)
        return sorted(set(gbas), key=lambda gba: roman_numerals.get(gba, 100))

    def __loop_classification(self) -> str:
        if not self.loops or len(self.loops) != 3 or any([loop.loop_type == 'n/a' for loop in self.loops]):
            return 'n/a'
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
        fingerprint = ''.join([loop.loop_type[0] for loop in self.loops])
        if fingerprint not in loop_classes:
            logging.error(f'Unknown loop classification: {fingerprint}')
            return 'n/a'
        subtype = 'a' if self.loops[0 if fingerprint != 'dpd' else 1].loop_type[-1] == '-' else 'b'
        return loop_classes[fingerprint] + subtype

    def __tracts(self) -> List[Tract]:
        if not self.tetrad_pairs:
            return [Tract([self.tetrads[0].nucleotides[i]]) for i in range(4)]

        tracts = [[self.tetrad_pairs[0].tetrad1.nucleotides[i]] for i in range(4)]
        for tetrad_pair in self.tetrad_pairs:
            for i in range(4):
                tracts[i].append(tetrad_pair.stacked[tracts[i][-1]])
        return [Tract(nts) for nts in tracts]

    def __loops(self) -> List[Loop]:
        if len(self.tetrads) == 1:
            return []

        loops = []
        tetrad_nucleotides = sorted([nt for tetrad in self.tetrads for nt in tetrad.nucleotides])
        for i in range(1, len(tetrad_nucleotides)):
            nprev = tetrad_nucleotides[i - 1]
            ncur = tetrad_nucleotides[i]
            if ncur.index - nprev.index > 1 and ncur.chain == nprev.chain:
                for tract in self.tracts:
                    if nprev in tract.nucleotides and ncur in tract.nucleotides:
                        break
                else:
                    nts = list(filter(lambda nt: nprev.index < nt.index < ncur.index, self.nucleotides.values()))
                    loop_type = self.__loop_type(nprev, ncur)
                    loops.append(Loop(nts, loop_type))
        return loops

    def __loop_type(self, nt1: Nucleotide, nt2: Nucleotide):

        def is_anticlockwise(n1: Nucleotide, n2: Nucleotide, t: Tetrad):
            for p1 in t.pairs:
                p2 = p1.reverse()
                if (n1, n2) == p1.pair:
                    return p1.score() < p2.score()
                if (n1, n2) == p2.pair:
                    return p1.score() > p2.score()
            logging.error(f'Failed to find a pair between {n1} and {n2}')
            return None

        for tetrad in self.tetrads:
            if {nt1, nt2}.issubset(tetrad.set):
                for pair in map(lambda p: {*p.pair}, tetrad.pairs):
                    if {nt1, nt2} == pair:
                        return 'lateral' + ('-' if is_anticlockwise(nt1, nt2, tetrad) else '+')
                return 'diagonal'

        tract_with_nt2 = next(filter(lambda tract: nt2 in tract.nucleotides, self.tracts))

        for nt2 in tract_with_nt2.nucleotides:
            pair = frozenset([nt1, nt2])
            for tetrad in self.tetrads:
                if nt2 in tetrad.nucleotides and any([frozenset(nt_pair.pair) == pair for nt_pair in tetrad.pairs]):
                    return 'propeller' + ('-' if is_anticlockwise(nt1, nt2, tetrad) else '+')

        logging.warning(f'Failed to classify the loop between {nt1} and {nt2}')
        return 'n/a'


class Helix:

    def __init__(self, tetrads: List[Tetrad], tetrad_pairs: List[TetradPair], nucleotides: Dict[str, Nucleotide]):
        self.tetrads: List[Tetrad] = tetrads
        self.tetrad_pairs: List[TetradPair] = tetrad_pairs
        self.nucleotides: Dict[str, Nucleotide] = nucleotides
        self.quadruplexes: List[Quadruplex] = self.__quadruplexes()

    def __iter__(self):
        return iter(self.tetrads)

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

    def __quadruplexes(self):
        if not self.tetrad_pairs:
            return [Quadruplex(self.tetrads, [], self.nucleotides)]

        quadruplexes = list()
        tetrads = list()
        for tetrad in [self.tetrad_pairs[0].tetrad1] + [tetrad_pair.tetrad2 for tetrad_pair in self.tetrad_pairs]:
            if tetrads:
                if tetrad.chains.isdisjoint(tetrads[-1].chains):
                    quadruplexes.append(Quadruplex(tetrads, self.__filter_tetrad_pairs(tetrads), self.nucleotides))
                    tetrads = list()
            tetrads.append(tetrad)

        quadruplexes.append(Quadruplex(tetrads, self.__filter_tetrad_pairs(tetrads), self.nucleotides))

        return quadruplexes

    def __filter_tetrad_pairs(self, tetrads: Iterable[Tetrad]) -> List[TetradPair]:
        chains = set()
        for tetrad in tetrads:
            chains.update(tetrad.chains)

        def check_tetrad(t: Tetrad) -> bool:
            return not t.chains.isdisjoint(chains)

        def check_pair(tp: TetradPair) -> bool:
            return check_tetrad(tp.tetrad1) and check_tetrad(tp.tetrad2)

        return list(filter(check_pair, self.tetrad_pairs))


class Analysis:

    def __init__(self, data: dict, structure3d: Optional[Structure3D]):
        self.nucleotides: Dict[str, Nucleotide] = {nt['nt_id']: Nucleotide(nt, structure3d) for nt in data['nts']}
        self.stacking: Set[Tuple[Nucleotide]] = self.__read_stacking(data)
        self.pairs: Dict[Tuple[Nucleotide, Nucleotide], Pair] = self.__read_pairs(data)
        self.canonical: Set[Pair] = set(filter(Pair.is_canonical, self.pairs.values()))
        self.metal_ions: List[Atom3D] = structure3d.find_metal_ions() if structure3d else []
        self.graph: Dict[Nucleotide, List[Nucleotide]] = {}
        self.tetrads: List[Tetrad] = []
        self.tetrad_pairs: List[TetradPair] = []
        self.stems: Dict[Tetrad, List[Tetrad]] = {}
        self.stacks: Dict[Tetrad, List[Tetrad]] = {}
        self.helices: List[Helix] = []
        self.sequence: str
        self.line1: str
        self.line2: str
        self.shifts: Dict[Nucleotide, int]

    def __str__(self):
        builder = f'Chain order: {", ".join(self.chain_order())}\n'
        for helix in self.helices:
            builder += str(helix)
        builder += f'{self.sequence}\n{self.line1}\n{self.line2}'
        return builder

    def build_graph(self, strict: bool):
        graph = defaultdict(list)
        for pair in self.pairs.values():
            if strict and pair.lw not in ('cWH', 'cHW'):
                continue
            nt1, nt2 = pair.pair
            graph[nt1].append(nt2)
        self.graph = graph

    def find_tetrads(self, no_reorder=False):
        # search for a tetrad: i -> j -> k -> l
        #                      ^--------------^
        tetrads = set()
        for i in self.graph:
            for j in filter(lambda x: x != i, self.graph[i]):
                for k in filter(lambda x: x not in (i, j), self.graph[j]):
                    for l in filter(lambda x: x not in (i, j, k) and x in self.graph[i], self.graph[k]):
                        if Tetrad.is_valid(i, j, k, l, self.pairs):
                            tetrads.add(Tetrad(i, j, k, l, self.pairs, no_reorder))

        # build graph of tetrads
        while True:
            graph = {t: [] for t in tetrads}
            for (ti, tj) in itertools.combinations(tetrads, 2):
                if not ti.is_disjoint(tj):
                    graph[ti].append(tj)
                    graph[tj].append(ti)

            # remove tetrad which conflicts the most with others
            # in case of a tie, remove one which has the worst planarity deviation
            candidates = sorted(tetrads, key=lambda t: (len(graph[t]), t.planarity_deviation), reverse=True)
            if len(graph[candidates[0]]) > 0:
                tetrads.remove(candidates[0])
            else:
                break

        self.tetrads = sorted(tetrads, key=lambda t: min(t.nucleotides[i].index for i in range(4)))

    def find_stacks(self, stacking_mismatch: int):
        stackings = defaultdict(list)
        for ti, tj in itertools.combinations(self.tetrads, 2):
            i, j = ti.count_non_stacked_bases(tj, self.stacking)
            if i <= stacking_mismatch and j <= stacking_mismatch:
                stackings[ti].append(tj)
                stackings[tj].append(ti)
        self.stacks = stackings

    def find_tetrad_pairs_and_helices(self, stacking_mismatch: int):

        def is_next_by_stacking(nt1, nt2):
            for stack in self.stacking:
                if nt1 in stack and nt2 in stack:
                    return abs(stack.index(nt2) - stack.index(nt1)) == 1
            return False

        def is_next_sequentially(nt1, nt2):
            return nt1.chain == nt2.chain and abs(nt1.index - nt2.index) == 1

        tetrad_scores = defaultdict(dict)

        for ti, tj in itertools.combinations(self.tetrads, 2):
            nts1 = ti.nucleotides
            best_score = 0
            best_order = tj.nucleotides

            for nts2 in itertools.permutations(tj.nucleotides):
                score = 0
                for i in range(4):
                    if is_next_by_stacking(nts1[i], nts2[i]) or is_next_sequentially(nts1[i], nts2[i]):
                        score += 1
                if score > best_score:
                    best_score = score
                    best_order = nts2
                if best_score == 4:
                    break

            tetrad_scores[ti][tj] = (best_score, nts1, best_order)
            tetrad_scores[tj][ti] = (best_score, best_order, nts1)

        tetrads = list(self.tetrads)
        best_score = 0
        best_order = tetrads

        for ti in tetrads:
            score = 0
            order = [ti]
            candidates = set(self.tetrads) - {ti}

            while candidates:
                tj = max([tj for tj in candidates], key=lambda tk: tetrad_scores[ti][tk][0])
                score += tetrad_scores[ti][tj][0]
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
            score = tetrad_scores[ti][tj][0]

            if score > (4 - stacking_mismatch):
                nts1, nts2 = tetrad_scores[ti][tj][1:]
                stacked = dict([(nts1[i], nts2[i]) for i in range(4)])
                stacked.update({v: k for k, v in stacked.items()})
                tetrad_pairs.append(TetradPair(ti, tj, stacked))

        helices = []
        helix_tetrads = []
        helix_tetrad_pairs = []

        for tp in tetrad_pairs:
            ti, tj = tp.tetrad1, tp.tetrad2
            if not helix_tetrads:
                helix_tetrads.append(ti)
            score = tetrad_scores[helix_tetrads[-1]][tj][0]
            if score >= (4 - stacking_mismatch):
                helix_tetrads.append(tj)
                helix_tetrad_pairs.append(tp)
            else:
                helices.append(Helix(helix_tetrads, helix_tetrad_pairs, self.nucleotides))
                helix_tetrads = [ti, tj]
                helix_tetrad_pairs = [tp]

        if helix_tetrads:
            helices.append(Helix(helix_tetrads, helix_tetrad_pairs, self.nucleotides))

        for tetrad in self.tetrads:
            if not any([tetrad in helix.tetrads for helix in helices]):
                helices.append(Helix([tetrad], [], self.nucleotides))

        self.tetrad_pairs = tetrad_pairs
        self.helices = helices

    def find_best_chain_reorder(self, stacking_mismatch: int):
        scores = {'O+': 0, 'O-': 1, 'N+': 2, 'N-': 3, 'Z+': 4, 'Z-': 5}
        chain_groups = self.__group_related_chains()
        final_order = []

        for chains in filter(lambda x: len(x) > 1, chain_groups):
            best_permutation, best_score = chains, (1e10, 1e10)
            for permutation in sorted(itertools.permutations(chains)):
                self.__reorder_chains(permutation)
                classification = list(self.__get_classification())
                score = (sum(scores[c] for c in classification), self.__sum_squares_chain_distances(permutation))
                logging.debug(f'Checking reorder: {" ".join(permutation)} {"".join(classification)}')
                if score < best_score:
                    best_score = score
                    best_permutation = permutation
            final_order.extend(best_permutation)

        if final_order:
            self.__reorder_chains(final_order)
            self.tetrads = sorted(self.tetrads, key=lambda t: min(nt.index for nt in t.nucleotides))
            self.find_tetrad_pairs_and_helices(stacking_mismatch)
            logging.debug(f'Selected reorder: {" ".join(final_order)} {"".join(list(self.__get_classification()))}')

    def chain_order(self) -> Dict[str, int]:
        return {nt.model_chain: 0 for nt in sorted(self.nucleotides.values(), key=lambda nt: nt.index)}

    def analyze_metal_ions(self):
        for ion in self.metal_ions:
            min_distance = float('inf')
            min_tetrad = None
            min_nt = None
            channel = None

            for tetrad in self.tetrads:
                distance = numpy.linalg.norm(ion.coordinates() - tetrad.center)
                if distance < min_distance:
                    min_distance = distance
                    min_tetrad = tetrad
                    min_nt = None
                    channel = True

                for nt in tetrad.nucleotides:
                    distance = numpy.linalg.norm(ion.coordinates() - nt.outermost_atom().coordinates())
                    if distance < min_distance:
                        min_distance = distance
                        min_tetrad = tetrad
                        min_nt = nt
                        channel = False

            # TODO: verify threshold of 6.0 Angstroms
            if min_distance < 6.0:
                if channel:
                    min_tetrad.ions_channel.append(ion)
                else:
                    if min_nt not in min_tetrad.ions_outside:
                        min_tetrad.ions_outside[min_nt] = []
                    min_tetrad.ions_outside[min_nt].append(ion)
            else:
                logging.debug(f'Skipping an ion, because it is too far from any tetrad (distance={min_distance})')

    def compute_twoline_dotbracket(self):
        layer1, layer2 = [], []
        for tetrad in self.tetrads:
            layer1.extend((tetrad.pairs[0], tetrad.pairs[2]))
            layer2.extend((tetrad.pairs[1], tetrad.pairs[3]))
        self.sequence, self.line1, self.shifts = self.__elimination_conflicts(layer1)
        self.line2 = self.__elimination_conflicts(layer2)[1]

    def __reorder_chains(self, chain_order: Iterable[str]):
        i = 1
        for chain in chain_order:
            for nt in self.nucleotides.values():
                if nt.model_chain == chain:
                    nt.index = i
                    i += 1
        for nt in self.nucleotides.values():
            if nt.model_chain not in chain_order:
                nt.index = i
                i += 1
        for tetrad in self.tetrads:
            tetrad.reorder()

    def __group_related_chains(self):
        chains = dict()
        for q in self.helices:
            chains[q] = {n.model_chain for t in q.tetrads for n in t.nucleotides}
        changed = True
        while changed:
            changed = False
            for q1, q2 in itertools.combinations(chains.keys(), 2):
                if not chains[q1].isdisjoint(chains[q2]):
                    chains[q1].update(chains[q2])
                    del chains[q2]
                    changed = True
                    break
        return chains.values()

    def __sum_squares_chain_distances(self, chain_order: Tuple) -> int:
        chain_pairs = set(
            (frozenset((p.pair[0].chain, p.pair[1].chain)) for h in self.helices for t in h.tetrads
             for p in t.pairs
             if p.pair[0].chain != p.pair[1].chain and p.pair[0].chain in chain_order and p.pair[1].chain in chain_order
            ))
        index = {chain: chain_order.index(chain) for chain in chain_order}
        sum_sq = 0
        for c1, c2 in chain_pairs:
            sum_sq += (index[c1] - index[c2])**2
        return sum_sq

    def __get_classification(self):
        for q in self.helices:
            for t in q.tetrads:
                yield t.get_classification()

    def __read_stacking(self, data: dict) -> Set[Tuple[Nucleotide]]:
        stacking = set()
        for stack in data['stacks']:
            stack = tuple(filter(lambda nt_id: nt_id in self.nucleotides, stack['nts_long'].split(',')))
            if len(stack) > 1:
                stacking.add(tuple(self.nucleotides[nt] for nt in stack))
        return stacking

    def __read_pairs(self, data: dict) -> Dict[Tuple[Nucleotide, Nucleotide], Pair]:
        pairs = dict()
        for pair in data['pairs']:
            nt1, nt2, lw, saenger = pair['nt1'], pair['nt2'], pair['LW'], pair['Saenger']
            if Pair.is_valid(lw):
                nt1, nt2 = self.nucleotides[nt1], self.nucleotides[nt2]
                pair = Pair(nt1, nt2, lw, saenger)
                pairs[(nt1, nt2)] = pair
                pairs[(nt2, nt1)] = pair.reverse()
        return pairs

    def __elimination_conflicts(self, pairs):
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
                pair, _ = sorted(conflicts.items(), key=lambda x: (len(x[1]), x[0].pair[0]), reverse=True)[0]
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
            nt1, nt2 = sorted(pair.pair)
            dotbracket[nt1] = opening[order]
            dotbracket[nt2] = closing[order]

        sequence = ''
        structure = ''
        shifts = dict()
        shift_value = 0
        chain = None
        for nt in sorted(self.nucleotides.values()):
            if chain and chain != nt.model_chain:
                sequence += '-'
                structure += '-'
                shift_value += 1
            sequence += nt.short_name
            structure += dotbracket.get(nt, '.')
            shifts[nt] = shift_value
            chain = nt.model_chain
        return sequence, structure, shifts


class Visualizer:

    def __init__(self,
                 tetrads: Iterable[Tetrad],
                 sequence: str,
                 shifts: Dict[Nucleotide, int],
                 nucleotides: Iterable[Nucleotide] = tuple(),
                 canonical: Iterable[Pair] = tuple()):
        self.tetrads = tetrads
        self.sequence = sequence
        self.shifts = shifts
        self.nucleotides = nucleotides if nucleotides else self.__extract_nucleotides()
        self.canonical = canonical
        self.onz = {pair: tetrad.get_classification() for tetrad in tetrads for pair in tetrad.pairs}

    def visualize(self, prefix: str, suffix: str):
        tempdir = os.path.join(tempfile.gettempdir(), '.tetrads')
        os.makedirs(tempdir, exist_ok=True)

        fd, fasta = tempfile.mkstemp('.fasta', dir=tempdir)
        os.close(fd)
        with open(fasta, 'w') as fastafile:
            fastafile.write(f'>{prefix}-{suffix}\n')
            fastafile.write(''.join(self.sequence))

        layer1, layer2 = [], []
        for tetrad in self.tetrads:
            layer1.extend((tetrad.pairs[0], tetrad.pairs[2]))
            layer2.extend((tetrad.pairs[1], tetrad.pairs[3]))
        helix1 = self.__to_helix(layer1, tempdir, self.canonical)
        helix2 = self.__to_helix(layer2, tempdir)

        currdir = os.path.dirname(os.path.realpath(__file__))
        output_pdf = f'{prefix}-{suffix}.pdf'
        run = subprocess.run([os.path.join(currdir, 'quadraw.R'), fasta, helix1, helix2, output_pdf],
                             stdout=subprocess.PIPE,
                             stderr=subprocess.PIPE)
        if run.returncode == 0:
            print('\nPlot:', output_pdf)
        else:
            logging.error(f'Failed to prepare visualization, reason:\n  {run.stderr.decode()}')

    def __extract_nucleotides(self):
        nucleotides = list()
        for tetrad in self.tetrads:
            nucleotides.extend(tetrad.nucleotides)
        return nucleotides

    def __to_helix(self, layer: List[Pair], tempdir: str, canonical: Iterable[Pair] = tuple()) -> str:
        fd, name = tempfile.mkstemp('.helix', dir=tempdir)
        os.close(fd)
        onz_value = {'O+': 1, 'O-': 2, 'N+': 3, 'N-': 4, 'Z+': 5, 'Z-': 6, 'n/a': 7}
        nucleotides = sorted(self.nucleotides)

        with open(name, 'w') as helixfile:
            helixfile.write(f'#{len(self.sequence) + 1}\n')
            helixfile.write('i\tj\tlength\tvalue\n')
            for pair in layer:
                x, y = pair.pair
                x, y = nucleotides.index(x) + 1 + self.shifts[x], nucleotides.index(y) + 1 + self.shifts[y]
                onz = self.onz[pair]
                helixfile.write(f'{x}\t{y}\t1\t{onz_value[onz]}\n')
            for pair in canonical:
                x, y = pair.pair
                x, y = nucleotides.index(x) + 1 + self.shifts[x], nucleotides.index(y) + 1 + self.shifts[y]
                helixfile.write(f'{x}\t{y}\t1\t8\n')
        return name


class AnalysisSimple:

    def __init__(self, data: dict):
        self.pairs = read_pairs(data)
        self.graph = dict()

    def build_graph(self):
        graph = defaultdict(list)
        for pair in self.pairs.values():
            nt1, nt2 = pair.pair
            graph[nt1].append(nt2)
        self.graph = graph

    def is_valid_tetrad(self, nt1, nt2, nt3, nt4):
        p1, p2, p3, p4 = self.pairs[(nt1, nt2)], self.pairs[(nt2, nt3)], self.pairs[(nt3, nt4)], self.pairs[(nt4, nt1)]
        for pi, pj in ((p1, p4), (p2, p1), (p3, p2), (p4, p3)):
            if pi.lw[1] == pj.lw[2]:
                return False
        return True

    def has_tetrads(self):
        tetrads = set()
        for i in self.graph:
            for j in filter(lambda x: x != i, self.graph[i]):
                for k in filter(lambda x: x not in (i, j), self.graph[j]):
                    for l in filter(lambda x: x not in (i, j, k) and x in self.graph[i], self.graph[k]):
                        if self.is_valid_tetrad(i, j, k, l):
                            tetrads.add(frozenset([i, j, k, l]))
                        if len(tetrads) > 1:
                            return True
        return False


def read_pairs(data: dict):
    pairs = dict()
    for pair in data['pairs']:
        nt1, nt2, lw, saenger = pair['nt1'], pair['nt2'], pair['LW'], pair['Saenger']
        if Pair.is_valid(lw):
            pair = Pair(nt1, nt2, lw, saenger)
            pairs[(nt1, nt2)] = pair
            pairs[(nt2, nt1)] = pair.reverse()
    return pairs


def center_of_mass(atoms):
    coords = [atom.coordinates() for atom in atoms]
    xs = (coord[0] for coord in coords)
    ys = (coord[1] for coord in coords)
    zs = (coord[2] for coord in coords)
    return numpy.array((sum(xs) / len(coords), sum(ys) / len(coords), sum(zs) / len(coords)))


def eltetrado(dssr: Dict, structure3d: Structure3D, strict: bool, no_reorder: bool, stacking_mismatch: int) -> Analysis:
    if 'pairs' not in dssr:
        return Analysis({}, None)

    structure = Analysis(dssr, structure3d)
    structure.build_graph(strict)
    structure.find_tetrads(no_reorder)
    structure.analyze_metal_ions()
    structure.find_stacks(stacking_mismatch)
    structure.find_tetrad_pairs_and_helices(stacking_mismatch)

    if not no_reorder:
        structure.find_best_chain_reorder(stacking_mismatch)

    structure.compute_twoline_dotbracket()

    return structure


def has_tetrad(dssr: Dict) -> bool:
    if 'pairs' not in dssr:
        return False

    structure = AnalysisSimple(dssr)
    structure.build_graph()
    return structure.has_tetrads()
