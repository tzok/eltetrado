import argparse
import gzip
import json
import os
import shutil
import subprocess
import sys
import tempfile

from collections import defaultdict

__version__ = '1.2.0'


class Pair:
    def __init__(self, nt1, nt2, lw):
        self.pair = (nt1, nt2)
        self.lw = lw

    def reverse(self):
        lw = '{}{}{}'.format(self.lw[0], self.lw[2], self.lw[1])
        return Pair(self.pair[1], self.pair[0], lw)


class Structure:
    def __init__(self, data: dict):
        self.pairs = self._read_pairs(data)
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

    def _read_pairs(self, data: dict):
        pairs = dict()
        for pair in data['pairs']:
            nt1, nt2 = pair['nt1'], pair['nt2']
            pair = Pair(nt1, nt2, pair['LW'])
            pairs[(nt1, nt2)] = pair
            pairs[(nt2, nt1)] = pair.reverse()
        return pairs


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--json', action='store_true', help='parse input file from JSON format')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(__version__))
    parser.add_argument('input',
                        help='a JSON file produced by DSSR if "--json" is used, otherwise a PDB or PDBx/MMCIF file to '
                             'be analyzed first by DSSR')
    return parser.parse_args()


def load_dssr_results(args):
    is_gzip = os.path.splitext(args.input)[1] == '.gz'

    if args.json:
        open_proxy = gzip.open if is_gzip else open
        with open_proxy(args.input) as jsonfile:
            dssr = jsonfile.read()
    else:
        tempdir = tempfile.mkdtemp()

        if is_gzip:
            fd, inputfile = tempfile.mkstemp(dir=tempdir)
            os.close(fd)
            with gzip.open(args.input) as compressed:
                with open(inputfile, 'wb') as uncompressed:
                    shutil.copyfileobj(compressed, uncompressed)
        else:
            inputfile = args.input

        currdir = os.path.dirname(os.path.realpath(__file__))
        shutil.copy(os.path.join(currdir, 'x3dna-dssr'), tempdir)
        dssr = subprocess.Popen(
            ['./x3dna-dssr', '-i={}'.format(os.path.abspath(inputfile)), '--json'],
            stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=tempdir)
        dssr, _ = dssr.communicate()

        shutil.rmtree(tempdir)

    try:
        return json.loads(dssr)
    except json.JSONDecodeError:
        print('Invalid JSON in', args.input, file=sys.stderr)
        exit(1)


def main():
    args = parse_arguments()
    dssr = load_dssr_results(args)

    if 'pairs' not in dssr:
        exit(1)

    structure = Structure(dssr)
    structure.build_graph()

    if structure.has_tetrads():
        exit()
    else:
        exit(1)

if __name__ == '__main__':
    main()
