import gzip
import logging
import os
import shutil
import subprocess
import tempfile

from typing import Dict, List, Optional

import orjson

from pdbx import PdbxReader

from eltetrado.compute.model import Atom3D, Structure3D, Residue3D


def run_dssr(pdb_path: str) -> Dict:
    # check if x3dna-dssr is on $PATH
    dssr_binary = shutil.which('x3dna-dssr')

    if dssr_binary is None:
        logging.error('Missing x3dna-dssr on $PATH, please install the application')
        return {}

    tempdir = tempfile.mkdtemp()

    try:
        shutil.copy(dssr_binary, tempdir)
        dssr = subprocess.Popen(['./x3dna-dssr', '-i={}'.format(os.path.abspath(pdb_path)), '--json'],
                                stdout=subprocess.PIPE,
                                stderr=subprocess.PIPE,
                                cwd=tempdir)
        dssr, _ = dssr.communicate()
        return orjson.loads(dssr.decode('utf-8'))
    finally:
        shutil.rmtree(tempdir)


def load_dssr_results(json_path: Optional[str] = None, pdb_path: Optional[str] = None) -> Dict:
    # load provided JSON file
    if json_path:
        with open(json_path) as jsonfile:
            return orjson.loads(jsonfile.read())

    # analyze PDB or PDBx/mmCIF file now
    if pdb_path:
        # if not gzipped, run directly on the file
        root, ext = os.path.splitext(pdb_path)

        if ext != '.gz':
            return run_dssr(pdb_path)

        # otherwise ungzip first
        fd, ungzipped = tempfile.mkstemp(suffix=os.path.splitext(root)[1])
        os.close(fd)

        try:
            with gzip.open(pdb_path, 'rb') as infile:
                with open(ungzipped, 'wb') as outfile:
                    outfile.write(infile.read())
            return run_dssr(ungzipped)
        finally:
            os.remove(ungzipped)

    logging.error('Neither DSSR JSON or PDB / mmCIF path supplied')
    return {}


def group_atoms(atoms: List[Atom3D]) -> Structure3D:
    if not atoms:
        return Structure3D([])

    key_previous = (atoms[0].chain_identifier, atoms[0].residue_number, atoms[0].insertion_code)
    residue_atoms = [atoms[0]]
    residues = []

    for atom in atoms[1:]:
        key = (atom.chain_identifier, atom.residue_number, atom.insertion_code)
        if key == key_previous:
            residue_atoms.append(atom)
        else:
            residues.append(
                Residue3D(residue_atoms, residue_atoms[0].residue_name, residue_atoms[0].chain_identifier,
                          residue_atoms[0].residue_number, residue_atoms[0].insertion_code))
            key_previous = key
            residue_atoms = [atom]

    residues.append(
        Residue3D(residue_atoms, residue_atoms[0].residue_name, residue_atoms[0].chain_identifier,
                  residue_atoms[0].residue_number, residue_atoms[0].insertion_code))
    return Structure3D(residues)


def parse_pdb(inputfile: str, is_gzip: bool) -> Structure3D:
    fd = gzip.open(inputfile, 'rt') if is_gzip else open(inputfile)
    atoms = []

    try:
        for line in fd:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                alternate_location = line[16]
                if alternate_location != ' ':
                    continue
                atom_name = line[12:16].strip()
                residue_name = line[18:20].strip()
                chain_identifier = line[21]
                residue_number = int(line[22:26].strip())
                insertion_code = line[26]
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                atoms.append(Atom3D(atom_name, residue_name, chain_identifier, residue_number, insertion_code, x, y, z))
    finally:
        fd.close()

    return group_atoms(atoms)


def parse_cif(inputfile: str, is_gzip: bool) -> Structure3D:
    fd = gzip.open(inputfile, 'rt') if is_gzip else open(inputfile)

    try:
        data = []
        reader = PdbxReader(fd)
        reader.read(data)
    finally:
        fd.close()

    atoms = []

    if data:
        atom_site = data[0].get_object('atom_site')

        for row in atom_site.row_list:
            row_dict = dict(zip(atom_site.attribute_list, row))
            atom_name = row_dict['auth_atom_id']
            residue_name = row_dict['auth_comp_id']
            chain_identifier = row_dict['auth_asym_id']
            residue_number = int(row_dict['auth_seq_id'])
            insertion_code = row_dict['pdbx_PDB_ins_code'] if row_dict['pdbx_PDB_ins_code'] else ' '
            x = float(row_dict['Cartn_x'])
            y = float(row_dict['Cartn_y'])
            z = float(row_dict['Cartn_z'])
            atoms.append(Atom3D(atom_name, residue_name, chain_identifier, residue_number, insertion_code, x, y, z))

    return group_atoms(atoms)


def read_3d_structure(inputfile: str) -> Optional[Structure3D]:
    if not inputfile:
        return None

    root, ext = os.path.splitext(inputfile)

    if ext == '.gz':
        is_gzip = True
        root, ext = os.path.splitext(root)
    else:
        is_gzip = False

    if ext == '.pdb':
        return parse_pdb(inputfile, is_gzip)
    elif ext == '.cif':
        return parse_cif(inputfile, is_gzip)
    else:
        logging.error(f'Unknown file type: {inputfile}')
        return None
