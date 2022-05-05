import http
import os
from typing import List, TextIO, Tuple, Dict, Union, Optional

import requests
from mmcif.io import IoAdapter

from eltetrado.model import Atom3D, Residue3D, ResidueLabel, ResidueAuth, Structure3D, Structure2D

RNAPOLIS_WS_URL = os.getenv('RNAPOLIS_WS_URL', 'https://rnapolis-ws.cs.put.poznan.pl/api')


def read_2d_structure(cif_or_pdb: TextIO, model: int) -> Structure2D:
    cif_or_pdb.seek(0)
    result = requests.post(f'{RNAPOLIS_WS_URL}/analyze/{model}', cif_or_pdb.read(),
                           headers={'Content-Type': 'text/plain'}, timeout=60 * 5)
    if result.status_code != http.HTTPStatus.OK:
        result.raise_for_status()
    return Structure2D(**result.json())


def read_3d_structure(cif_or_pdb: TextIO, model: int) -> Structure3D:
    atoms, modified, sequence = parse_cif(cif_or_pdb) if is_cif(cif_or_pdb) else parse_pdb(cif_or_pdb)
    atoms = list(filter(lambda atom: atom.model == model, atoms))
    return group_atoms(atoms, modified, sequence)


def is_cif(cif_or_pdb: TextIO) -> bool:
    cif_or_pdb.seek(0)
    for line in cif_or_pdb.readlines():
        if line.startswith('_atom_site'):
            return True
    return False


def parse_cif(cif: TextIO) \
        -> Tuple[List[Atom3D], Dict[Union[ResidueLabel, ResidueAuth], str], Dict[Tuple[str, int], str]]:
    cif.seek(0)

    io_adapter = IoAdapter()
    data = io_adapter.readFile(cif.name)
    atoms = []
    modified = {}
    sequence = {}

    if data:
        atom_site = data[0].getObj('atom_site')
        mod_residue = data[0].getObj('pdbx_struct_mod_residue')
        entity_poly = data[0].getObj('entity_poly')

        if atom_site:
            for row in atom_site.getRowList():
                row_dict = dict(zip(atom_site.getAttributeList(), row))

                label_chain_name = row_dict.get('label_asym_id', None)
                label_residue_number = try_parse_int(row_dict.get('label_seq_id', None))
                label_residue_name = row_dict.get('label_comp_id', None)
                auth_chain_name = row_dict.get('auth_asym_id', None)
                auth_residue_number = try_parse_int(row_dict.get('auth_seq_id', None))
                auth_residue_name = row_dict.get('auth_comp_id', None)
                insertion_code = row_dict.get('pdbx_PDB_ins_code', None)

                if label_chain_name is None and auth_chain_name is None:
                    raise RuntimeError(f'Cannot parse an atom line with empty chain name: {row}')
                if label_residue_number is None and auth_residue_number is None:
                    raise RuntimeError(f'Cannot parse an atom line with empty residue number: {row}')
                if label_residue_name is None and auth_residue_name is None:
                    raise RuntimeError(f'Cannot parse an atom line with empty residue name: {row}')

                label = None
                if label_chain_name and label_residue_number and label_residue_name:
                    label = ResidueLabel(label_chain_name, label_residue_number, label_residue_name)

                auth = None
                if auth_chain_name and auth_residue_number and auth_residue_name and insertion_code:
                    auth = ResidueAuth(auth_chain_name, auth_residue_number, insertion_code, auth_residue_name)

                model = int(row_dict.get('pdbx_PDB_model_num', '1'))
                atom_name = row_dict['label_atom_id']
                x = float(row_dict['Cartn_x'])
                y = float(row_dict['Cartn_y'])
                z = float(row_dict['Cartn_z'])
                atoms.append(Atom3D(label, auth, model, atom_name, x, y, z))

        if mod_residue:
            for row in mod_residue.getRowList():
                row_dict = dict(zip(mod_residue.getAttributeList(), row))

                label_chain_name = row_dict.get('label_asym_id', None)
                label_residue_number = try_parse_int(row_dict.get('label_seq_id', None))
                label_residue_name = row_dict.get('label_comp_id', None)
                auth_chain_name = row_dict.get('auth_asym_id', None)
                auth_residue_number = try_parse_int(row_dict.get('auth_seq_id', None))
                auth_residue_name = row_dict.get('auth_comp_id', None)
                insertion_code = row_dict.get('PDB_ins_code', None)

                label = None
                if label_chain_name and label_residue_number and label_residue_name:
                    label = ResidueLabel(label_chain_name, label_residue_number, label_residue_name)

                auth = None
                if auth_chain_name and auth_residue_number and auth_residue_name and insertion_code:
                    auth = ResidueAuth(auth_chain_name, auth_residue_number, insertion_code, auth_residue_name)

                # TODO: is processing this data for each model separately required?
                # model = row_dict.get('PDB_model_num', '1')
                standard_residue_name = row_dict.get('parent_comp_id', 'n')

                modified[label] = standard_residue_name
                modified[auth] = standard_residue_name

        if entity_poly:
            for row in entity_poly.getRowList():
                row_dict = dict(zip(entity_poly.getAttributeList(), row))

                pdbx_strand_id = row_dict.get('pdbx_strand_id', None)
                pdbx_seq_one_letter_code_can = row_dict.get('pdbx_seq_one_letter_code_can', None)

                if pdbx_strand_id and pdbx_seq_one_letter_code_can:
                    for strand in pdbx_strand_id.split(','):
                        for i, letter in enumerate(pdbx_seq_one_letter_code_can):
                            sequence[(strand, i + 1)] = letter

    return atoms, modified, sequence


def parse_pdb(pdb: TextIO) -> Tuple[List[Atom3D], Dict[ResidueAuth, str], Dict]:
    pdb.seek(0)
    atoms = []
    modified = {}
    model = 1

    for line in pdb.readlines():
        if line.startswith('MODEL'):
            model = int(line[10:14].strip())
        elif line.startswith('ATOM') or line.startswith('HETATM'):
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
            auth = ResidueAuth(chain_identifier, residue_number, insertion_code, residue_name)
            atoms.append(Atom3D(None, auth, model, atom_name, x, y, z))
        elif line.startswith('MODRES'):
            original_name = line[12:15]
            chain_identifier = line[16]
            residue_number = int(line[18:22].strip())
            insertion_code = line[23]
            standard_residue_name = line[24:27].strip()
            auth = ResidueAuth(chain_identifier, residue_number, insertion_code, original_name)
            modified[auth] = standard_residue_name

    return atoms, modified, {}


def group_atoms(atoms: List[Atom3D], modified: Dict[Union[ResidueLabel, ResidueAuth], str],
                sequence: Dict[Tuple[str, int], str]) -> Structure3D:
    if not atoms:
        return Structure3D([])

    key_previous = (atoms[0].label, atoms[0].auth, atoms[0].model)
    residue_atoms = [atoms[0]]
    residues = []
    index = 1

    for atom in atoms[1:]:
        key = (atom.label, atom.auth, atom.model)
        if key == key_previous:
            residue_atoms.append(atom)
        else:
            label = key_previous[0]
            auth = key_previous[1]
            model = key_previous[2]
            name = get_residue_name(auth, label, modified)
            one_letter_name = get_one_letter_name(label, sequence, name)
            residues.append(Residue3D(index, name, one_letter_name, model, label, auth, tuple(residue_atoms)))
            index += 1
            key_previous = key
            residue_atoms = [atom]

    label = key_previous[0]
    auth = key_previous[1]
    model = key_previous[2]
    name = get_residue_name(auth, label, modified)
    one_letter_name = get_one_letter_name(label, sequence, name)
    residues.append(Residue3D(index, name, one_letter_name, model, label, auth, tuple(residue_atoms)))
    return Structure3D(residues)


def get_residue_name(auth: ResidueAuth, label: ResidueLabel,
                     modified: Dict[Union[ResidueAuth, ResidueLabel], str]) -> str:
    if auth in modified:
        name = modified[auth].lower()
    elif label in modified:
        name = modified[label].lower()
    elif auth:
        name = auth.name
    elif label:
        name = label.name
    else:
        # any nucleotide
        name = 'n'
    return name


def get_one_letter_name(label: ResidueLabel, sequence: Dict[Tuple[str, int], str], name: str) -> str:
    # try getting the value from _entity_poly first
    if label:
        key = (label.chain, label.number)
        if key in sequence:
            return sequence[key]

    # RNA
    if len(name) == 1:
        return name
    # DNA
    if len(name) == 2 and name[0].upper() == 'D':
        return name[1]
    # try the last letter of the name
    if str.isalpha(name[-1]):
        return name[-1]
    # any nucleotide
    return 'n'


def try_parse_int(s: str) -> Optional[int]:
    try:
        return int(s)
    except:
        return None
