import os
import itplz.mol_def
from collections import defaultdict


class Includes:
    def __init__(self):
        self.mapping = defaultdict(set)

    def add(self, file, reason):
        self.mapping[file].add(reason)

    def print(self):
        for f in self.mapping:
            print(f'#includes "{f}"; {self.mapping[f]}')


def create_includes(top, defs, molecules, out_file):
    """Given a list of molecules, generate the include statement for them."""

    incs = Includes()
    for mol in molecules:
        for idx, source in enumerate([top, defs]):
            if source.mol_defs.get(mol):
                for atom in source.mol_defs.get(mol).get_atom_names():
                    if source.atom_defs.get(atom):
                        incs.add(source.atom_defs.get(atom).file, "atomtypes")
                for atom in source.mol_defs.get(mol).get_bond_names():
                    if source.bond_types.get(atom):
                        incs.add(source.bond_types.get(atom).file, "bondtypes")
                for atom in source.mol_defs.get(mol).get_angle_names():
                    if source.angle_types.get(atom):
                        incs.add(source.angle_types.get(atom).file, "angletypes")

    for atom in top.get_atom_names():
        if defs.atom_defs.get(atom):
            incs.add(defs.atom_defs.get(atom).file, "atomtypes")

    for bond in top.get_bond_names():
        b = defs.bond_types.get(bond)
        if b is not None:
            incs.add(b.file, "bondtypes")

    for angle in top.get_angle_names():
        b = defs.angle_types.get(angle)
        if b is not None:
            incs.add(b.file, "angletypes")

    for mol in molecules:
        if top.mol_defs.get(mol) is not None:
            incs.add(top.mol_defs.get(mol).file, mol)
        if defs.mol_defs.get(mol) is not None:
            incs.add(defs.mol_defs.get(mol).file, mol)

    # for inc in incs.mapping:
    #    writer.write_include(inc, ",".join(incs.mapping[inc]))
    return incs


def run(top_file, topology, mdls, out_file, verbose=False):
    try:
        last_include = max(topology.include_lines.values())
    except Exception as e:
        last_include = 0
    incs = create_includes(topology, mdls, topology.get_molecule_names(), out_file)

    with open(out_file, "w") as f:
        for inc in incs.mapping:
            i_path = os.path.relpath(inc, os.path.dirname(out_file))
            f.write(f'#include "{i_path}" ; {", ".join(incs.mapping[inc])}\n')

        with open(top_file, "r") as top_f:
            for line in top_f.readlines():
                if not line.startswith("#include"):
                    f.write(line)
