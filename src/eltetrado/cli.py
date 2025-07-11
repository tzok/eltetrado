import argparse
import gzip
import logging
import os
import sys
import tempfile
from typing import IO, List, Optional

import orjson
import rnapolis.annotator
import rnapolis.parser
from rnapolis.adapter import ExternalTool, parse_external_output, auto_detect_tool
from rnapolis.annotator import BaseInteractions, LeontisWesthof
from rnapolis.common import BasePair, Residue, Stacking
from rnapolis.tertiary import Structure3D

from eltetrado.analysis import Visualizer, eltetrado, has_tetrad
from eltetrado.dto import generate_dto


def eltetrado_cli(args=sys.argv[1:]):
    with open(os.path.join(os.path.dirname(__file__), "VERSION")) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path to input PDB or PDBx/mmCIF file")
    parser.add_argument("-o", "--output", help="(optional) path for output JSON file")
    parser.add_argument(
        "-m", "--model", help="(optional) model number to process", default=1, type=int
    )
    parser.add_argument(
        "--no-reorder",
        action="store_true",
        help="chains of bi- and tetramolecular quadruplexes should be reordered to be able to have "
        "them classified; when this is set, chains will be processed in original order, which for "
        "bi-/tetramolecular means that they will likely be misclassified; use with care!",
    )
    parser.add_argument(
        "--complete-2d",
        action="store_true",
        help="when set, the visualization will also show canonical base pairs to provide context for "
        "the quadruplex",
    )
    parser.add_argument(
        "--image",
        metavar="DIR",
        help="directory where visualization files (PDF) will be saved; "
        "if omitted, no images are generated",
        default=None,
    )
    parser.add_argument(
        "-e",
        "--external-files",
        nargs="*",
        default=[],
        help="path(s) to external tool output file(s); if omitted ElTetrado will compute interactions itself",
    )
    parser.add_argument(
        "--tool",
        choices=[t.value for t in ExternalTool],
        help="name of the external tool that produced the files (auto-detected when not provided)",
    )
    parser.add_argument(
        "-v", "--version", action="version", version="%(prog)s {}".format(version)
    )
    args = parser.parse_args(args)

    if not args.input:
        print(parser.print_help())
        sys.exit(1)

    cif_or_pdb = handle_input_file(args.input)
    structure3d = rnapolis.parser.read_3d_structure(
        cif_or_pdb, args.model, nucleic_acid_only=False
    )

    external_files: List[str] = args.external_files
    selected_tool: Optional[ExternalTool] = (
        ExternalTool(args.tool) if args.tool else None
    )

    if selected_tool is None and external_files:
        selected_tool = auto_detect_tool(external_files)
        logging.info(f"Auto-detected external tool: {selected_tool.value}")

    if selected_tool:
        if not external_files and selected_tool == ExternalTool.MAXIT:
            external_files = [args.input]
        base_interactions = parse_external_output(
            external_files, selected_tool, structure3d
        )
    else:
        base_interactions = rnapolis.annotator.extract_base_interactions(
            structure3d, args.model
        )

    analysis = eltetrado(
        base_interactions,
        structure3d,
        args.no_reorder,
    )
    print(analysis)

    if args.image is not None:
        visualizer = Visualizer(
            analysis, analysis.tetrads, args.complete_2d, analysis.global_index
        )

        basename = os.path.basename(args.input)
        root, ext = os.path.splitext(basename)
        if ext == ".gz":
            root, ext = os.path.splitext(root)
        prefix = root
        suffix = "str"
        visualizer.visualize(prefix, suffix, args.image)

        for i, helix in enumerate(analysis.helices):
            hv = Visualizer(
                analysis, helix.tetrads, args.complete_2d, analysis.global_index
            )
            suffix = "h{}".format(i + 1)
            hv.visualize(prefix, suffix, args.image)

            for j, quadruplex in enumerate(helix.quadruplexes):
                qv = Visualizer(
                    analysis,
                    quadruplex.tetrads,
                    args.complete_2d,
                    analysis.global_index,
                )
                qv.visualize(prefix, "{}-q{}".format(suffix, j + 1), args.image)

                for k, tetrad in enumerate(quadruplex.tetrads):
                    tv = Visualizer(
                        analysis, [tetrad], args.complete_2d, analysis.global_index
                    )
                    tv.visualize(
                        prefix, "{}-q{}-t{}".format(suffix, j + 1, k + 1), args.image
                    )

    if args.output:
        dto = generate_dto(analysis)

        with open(args.output, "wb") as jsonfile:
            jsonfile.write(orjson.dumps(dto))


def has_tetrad_cli(args=sys.argv[1:]):
    with open(os.path.join(os.path.dirname(__file__), "VERSION")) as f:
        version = f.read().strip()

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="path to input PDB or PDBx/mmCIF file")
    parser.add_argument(
        "-m", "--model", help="(optional) model number to process", default=1, type=int
    )
    parser.add_argument(
        "-e",
        "--external-files",
        nargs="*",
        default=[],
        help="path(s) to external tool output file(s); if omitted ElTetrado will compute interactions itself",
    )
    parser.add_argument(
        "--tool",
        choices=[t.value for t in ExternalTool],
        help="name of the external tool that produced the files (auto-detected when not provided)",
    )
    parser.add_argument(
        "--version", action="version", version="%(prog)s {}".format(version)
    )
    args = parser.parse_args(args)

    if not args.input:
        print(parser.print_help())
        sys.exit(1)

    cif_or_pdb = handle_input_file(args.input)
    structure3d = rnapolis.parser.read_3d_structure(
        cif_or_pdb, args.model, nucleic_acid_only=False
    )

    external_files: List[str] = args.external_files
    selected_tool: Optional[ExternalTool] = (
        ExternalTool(args.tool) if args.tool else None
    )

    if selected_tool is None and external_files:
        selected_tool = auto_detect_tool(external_files)
        logging.info(f"Auto-detected external tool: {selected_tool.value}")

    if selected_tool:
        if not external_files and selected_tool == ExternalTool.MAXIT:
            external_files = [args.input]
        base_interactions = parse_external_output(
            external_files, selected_tool, structure3d
        )
    else:
        base_interactions = rnapolis.annotator.extract_base_interactions(
            structure3d, args.model
        )

    print(has_tetrad(base_interactions, structure3d))


def handle_input_file(path) -> IO[str]:
    root, ext = os.path.splitext(path)

    if ext == ".gz":
        root, ext = os.path.splitext(root)
        file = tempfile.NamedTemporaryFile("w+", suffix=ext)
        with gzip.open(path, "rt") as f:
            file.write(f.read())
            file.seek(0)
    else:
        file = tempfile.NamedTemporaryFile("w+", suffix=ext)
        with open(path) as f:
            file.write(f.read())
            file.seek(0)
    return file


def read_secondary_structure_from_dssr(
    structure3d: Structure3D, model: int, dssr_json_path: str
) -> BaseInteractions:
    base_pairs: List[BasePair] = []
    stackings: List[Stacking] = []

    with open(dssr_json_path) as f:
        dssr = orjson.loads(f.read())

    for result in dssr.get("models", []):
        if result.get("model", None) == model:
            dssr = result.get("parameters", {})
            break

    for pair in dssr.get("pairs", []):
        nt1 = match_dssr_name_to_residue(structure3d, pair.get("nt1", None))
        nt2 = match_dssr_name_to_residue(structure3d, pair.get("nt2", None))
        lw = match_dssr_lw(pair.get("LW", None))

        if nt1 is not None and nt2 is not None and lw is not None:
            base_pairs.append(BasePair(nt1, nt2, lw, None))

    for stack in dssr.get("stacks", []):
        nts = [
            match_dssr_name_to_residue(structure3d, nt)
            for nt in stack.get("nts_long", "").split(",")
        ]
        for i in range(1, len(nts)):
            nt1 = nts[i - 1]
            nt2 = nts[i]
            if nt1 is not None and nt2 is not None:
                stackings.append(Stacking(nt1, nt2, None))

    return BaseInteractions(base_pairs, stackings, [], [], [])


def match_dssr_name_to_residue(
    structure3d: Structure3D, nt_id: Optional[str]
) -> Optional[Residue]:
    if nt_id is not None:
        nt_id = nt_id.split(":")[-1]
        for residue in structure3d.residues:
            if residue.full_name == nt_id:
                return residue
        logging.warn(f"Failed to find residue {nt_id}")
    return None


def match_dssr_lw(lw: Optional[str]) -> Optional[LeontisWesthof]:
    return LeontisWesthof[lw] if lw in dir(LeontisWesthof) else None


if __name__ == "__main__":
    eltetrado_cli()
