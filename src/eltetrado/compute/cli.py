import argparse
import os
import sys

import orjson

from eltetrado.compute.analysis import Visualizer, eltetrado, has_tetrad
from eltetrado.compute.io import load_dssr_results, read_3d_structure
from eltetrado.compute.model import generate_dto


def eltetrado_cli():
    with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', help='path to input PDB or PDBx/mmCIF file')
    parser.add_argument('--dssr-json', help='path to input JSON file generated with `x3dna-dssr --json`')
    parser.add_argument('--output', help='(optional) path for output JSON file')
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
                        help='chains of bi- and tetramolecular quadruplexes are reordered to be able to have them '
                        'classified; when this is set, chains will be processed in original order and '
                        'bi-/tetramolecular quadruplexes will not be classified')
    parser.add_argument('--complete-2d',
                        action='store_true',
                        help='when set, the visualization will also show canonical base pairs to provide context for '
                        'the quadruplex')
    parser.add_argument('--no-image',
                        action='store_true',
                        help='when set, the visualization will not be created at all')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version))
    args = parser.parse_args()

    if not args.pdb and not args.dssr_json:
        print(parser.print_help())
        sys.exit(1)

    dssr = load_dssr_results(args.dssr_json, args.pdb)
    structure3d = read_3d_structure(args.pdb)
    structure = eltetrado(dssr, structure3d, args.strict, args.no_reorder, args.stacking_mismatch)
    print(structure)

    visualizer = Visualizer(structure.tetrads, structure.sequence, structure.shifts, structure.nucleotides.values(),
                            structure.canonical if args.complete_2d else tuple())

    if not args.no_image:
        inputname = args.pdb if args.pdb else args.dssr_json
        basename = os.path.basename(inputname)
        root, ext = os.path.splitext(basename)
        if ext == '.gz':
            root, ext = os.path.splitext(root)
        prefix = root
        suffix = 'str'
        visualizer.visualize(prefix, suffix)

        for i, helix in enumerate(structure.helices):
            hv = Visualizer(helix.tetrads, structure.sequence, structure.shifts)
            suffix = 'h{}'.format(i + 1)
            hv.visualize(prefix, suffix)

            for j, quadruplex in enumerate(helix.quadruplexes):
                qv = Visualizer(quadruplex.tetrads, structure.sequence, structure.shifts)
                qv.visualize(prefix, '{}-q{}'.format(suffix, j + 1))

                for k, tetrad in enumerate(quadruplex.tetrads):
                    tv = Visualizer([tetrad], structure.sequence, structure.shifts)
                    tv.visualize(prefix, '{}-q{}-t{}'.format(suffix, j + 1, k + 1))

    if args.output:
        dto = generate_dto(structure)

        with open(args.output, 'wb') as jsonfile:
            jsonfile.write(orjson.dumps(dto))


def has_tetrad_cli():
    with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument('--pdb', help='path to input PDB or PDBx/mmCIF file')
    parser.add_argument('--dssr-json', help='path to input JSON file generated with `x3dna-dssr --json`')
    parser.add_argument('--version', action='version', version='%(prog)s {}'.format(version))
    args = parser.parse_args()

    if not args.pdb and not args.dssr_json:
        print(parser.print_help())
        sys.exit(1)

    dssr = load_dssr_results(args.dssr_json, args.pdb)
    sys.exit(0 if has_tetrad(dssr) else 1)
