import sys
import os
import glob
from collections import defaultdict


def is_header(line, name):
    return line[0] == '[' and name in line


class MolDef:
    """
    A basic mapping of molecule definition to filename.
    """
    __slots__ = ('name', 'file', 'line_num')

    def __init__(self, name: str, file: str, line_num: int):
        self.name = name
        self.file = file
        self.line_num = line_num

class MolDefLoader:
    """
    Read a directory of ITP files and gather up all molecule definitions
    """
    def __init__(self):
        self.mol_defs = defaultdict(list)
        self.include_tree = defaultdict(list)
        self.include_lines = {}
        self.failed_includes = []

    def load_directory(self, dirname):
        for filename in glob.glob(f"{dirname}/*"):
            if not os.path.isfile(filename):
                continue
            try:
                self.load_itp_file(filename)
            except Exception as e:
                continue

    def check_uniqueness(self, molecules=[]):
        for mol in molecules:
            defs = self.mol_defs[mol]
            if len(defs) > 1:
                print(f"{mol} is declared more than once:")
                for d in defs:
                    s = sorted(d.atom_names)
                    print(f"   {d.file}:{d.line_num} {s}")
            if len(defs) == 0:
                print(f"{mol} is not declared!")

    def contains(self, molecule):
        return self.mol_defs.get(molecule) is not None

    def get(self, molecule):
        includes = self.mol_defs.get(molecule, [])
        if len(includes) == 0:
            return None
        return includes[0]

    def load_itp_file(self, filename, parent=None):
        curr_molecule_name = None
        curr_molecule_line_num = None
        in_header = False
        line_num = 0

        with open(filename, 'r') as f:
            for line in f.readlines():
                line_num += 1
                line = line.strip()
                skip = False
                if len(line) > 0:
                    if line[0] == '[':
                        if curr_molecule_name is not None:
                            self.mol_defs[curr_molecule_name].append(MolDef(curr_molecule_name, filename, curr_molecule_line_num)) 
                            curr_molecule_name = None
                            curr_molecule_line_num = None
                        in_header = False
                        skip = True
                    if is_header(line, 'moleculetype'):
                        curr_molecule_line_num = line_num
                        in_header = True
                    if skip:
                        continue
                    if in_header and line[0] != ';':
                        curr_molecule_name = line.strip().split()[0]
            # End of file, add any molecule we've been reading.
            if curr_molecule_name is not None:
                self.mol_defs[curr_molecule_name].append(MolDef(curr_molecule_name, filename, curr_molecule_line_num)) 
                curr_molecule_name = None
                curr_molecule_line_num = None
        self.load_includes(filename)

    def load_includes(self, filename):
        with open(filename, 'r')  as f:
            lines = f.readlines()
            for  idx, line in enumerate(lines):
                if line.startswith("#include"):
                    start, include_filename, end = line.split('"')
                    try:
                        self.include_lines[filename] = idx
                        self.include_tree[filename].append(include_filename)
                        self.load_itp_file(include_filename)
                    except Exception as e:
                        self.failed_includes.append(include_filename)


#mdl = MolDefLoader("resources/martini/*.itp")
#mdl.check_uniqueness(['POPE', 'POPI', 'DPPC', 'ARG'])
#mdl.create_includes(['POPE', 'POPI','DPPC', 'ARG'])

