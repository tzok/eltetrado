import argparse
import gzip
import os
import sys
import tempfile
from typing import TextIO

import orjson

from eltetrado.analysis import Visualizer, eltetrado, has_tetrad
from eltetrado.model import generate_dto
from eltetrado.structure import read_3d_structure, read_2d_structure


def eltetrado_cli():
    with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='path to input PDB or PDBx/mmCIF file')
    parser.add_argument('-o', '--output', help='(optional) path for output JSON file')
    parser.add_argument('-m', '--model', help='(optional) model number to process', default=1, type=int)
    parser.add_argument('--stacking-mismatch',
                        help='a perfect tetrad stacking covers 4 nucleotides; this option can be used with value 1 or '
                             '2 to allow this number of nucleotides to be non-stacked with otherwise well aligned '
                             'tetrad [default=2]',
                        default=2,
                        type=int)
    parser.add_argument('--strict',
                        action='store_true',
                        help='nucleotides in tetrad are found when linked only by cWH pairing')
    parser.add_argument('--no-reorder',
                        action='store_true',
                        help='chains of bi- and tetramolecular quadruplexes should be reordered to be able to have '
                             'them classified; when this is set, chains will be processed in original order, which for '
                             'bi-/tetramolecular means that they will likely be misclassified; use with care!')
    parser.add_argument('--complete-2d',
                        action='store_true',
                        help='when set, the visualization will also show canonical base pairs to provide context for '
                             'the quadruplex')
    parser.add_argument('--no-image',
                        action='store_true',
                        help='when set, the visualization will not be created at all')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version))
    args = parser.parse_args()

    if not args.input:
        print(parser.print_help())
        sys.exit(1)

    cif_or_pdb = handle_input_file(args.input)
    structure3d = read_3d_structure(cif_or_pdb, args.model)
    structure2d = read_2d_structure(cif_or_pdb, args.model)

    analysis = eltetrado(structure2d, structure3d, args.strict, args.no_reorder, args.stacking_mismatch)
    print(analysis)

    if not args.no_image:
        visualizer = Visualizer(analysis, analysis.tetrads, args.complete_2d)

        basename = os.path.basename(args.input)
        root, ext = os.path.splitext(basename)
        if ext == '.gz':
            root, ext = os.path.splitext(root)
        prefix = root
        suffix = 'str'
        visualizer.visualize(prefix, suffix)

        for i, helix in enumerate(analysis.helices):
            hv = Visualizer(analysis, helix.tetrads, args.complete_2d)
            suffix = 'h{}'.format(i + 1)
            hv.visualize(prefix, suffix)

            for j, quadruplex in enumerate(helix.quadruplexes):
                qv = Visualizer(analysis, quadruplex.tetrads, args.complete_2d)
                qv.visualize(prefix, '{}-q{}'.format(suffix, j + 1))

                for k, tetrad in enumerate(quadruplex.tetrads):
                    tv = Visualizer(analysis, [tetrad], args.complete_2d)
                    tv.visualize(prefix, '{}-q{}-t{}'.format(suffix, j + 1, k + 1))

    if args.output:
        dto = generate_dto(analysis)

        with open(args.output, 'wb') as jsonfile:
            jsonfile.write(orjson.dumps(dto))


def has_tetrad_cli():
    with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='path to input PDB or PDBx/mmCIF file')
    parser.add_argument('-m', '--model', help='(optional) model number to process', default=1, type=int)
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version))
    args = parser.parse_args()

    if not args.input:
        print(parser.print_help())
        sys.exit(1)

    cif_or_pdb = handle_input_file(args.input)
    structure2d = read_2d_structure(cif_or_pdb, args.model)
    structure3d = read_3d_structure(cif_or_pdb, args.model)
    flag = has_tetrad(structure2d, structure3d)
    sys.exit(0 if flag else 1)


def handle_input_file(path) -> TextIO:
    root, ext = os.path.splitext(path)

    if ext == '.gz':
        root, ext = os.path.splitext(root)
        file = tempfile.NamedTemporaryFile('w+', suffix=ext)
        with gzip.open(path, 'rt') as f:
            file.write(f.read())
            file.seek(0)
    else:
        file = tempfile.NamedTemporaryFile('w+', suffix=ext)
        with open(path) as f:
            file.write(f.read())
            file.seek(0)
    return file


if __name__ == '__main__':
    eltetrado_cli()
