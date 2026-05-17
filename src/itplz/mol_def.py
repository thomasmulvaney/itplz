from collections import defaultdict


class Loc:
    """A location of a definition"""

    __slots__ = ("file", "line_num", "source")


class Definition:
    """
    Instance of a definition of something.
    """

    __slots__ = ("name", "file", "line_num", "source")

    def __init__(self, name: str, file: str, line_num: int, source: str):
        self.name = name
        self.file = file
        self.line_num = line_num
        self.source = source

    def key(self):
        return self.name


class MolDef(Definition):
    def __init__(self, name: str, file: str, line_num: int, source: str):
        self.name = name
        self.atoms = Index()
        self.bonds = Index()
        self.angles = Index()
        self.file = file
        self.line_num = line_num
        self.source = source

    def get_atom_names(self):
        return [a.name for aas in self.atoms.values() for a in aas]

    def get_bond_names(self):
        return [a.name for aas in self.bonds.values() for a in aas]

    def get_angle_names(self):
        return [a.name for aas in self.angles.values() for a in aas]

    def charge(self):
        total = 0
        for atom in self.atoms.values():
            total += atom[0].charge
        return total


class AtomDef(Definition):
    pass


class BondType(Definition):
    pass


class AngleType(Definition):
    pass


class BondRef(Definition):
    """A bond that depends on a BondType definition.

    Looks like:

       int int string

    Where string is the reference to a bond type.
    """

    pass


class AngleRef(Definition):
    """A bond that depends on a BondType definition.

    Looks like:

       int int int string

    Where string is the reference to a angle type.
    """

    pass


class Mol(Definition):
    __slots__ = ("molecule", "count", "file", "line_num", "source")

    def __init__(self, name: str, count: int, file: str, line_num: int, source: str):
        self.name = name
        self.count = count
        self.file = file
        self.line_num = line_num
        self.source = source


class Atom(Definition):
    __slots__ = ("id", "name", "charge", "file", "line_num", "source")

    def __init__(
        self, id: int, name: str, charge: float, file: str, line_num: int, source: str
    ):
        self.id = id
        self.name = name
        self.charge = charge
        self.file = file
        self.line_num = line_num
        self.source = source

    def key(self):
        return (self.id, self.name)


class Index:
    """Store sets of types for lookup."""

    def __init__(self):
        self._index = defaultdict(list)

    def add(self, item):
        self._index[item.key()].append(item)
        return item

    def keys(self):
        return self._index.keys()

    def values(self):
        return self._index.values()

    def get_all(self, item_name, source=None):
        return self._index.get(item_name)

    def get(self, item_name, source=None):
        vals = self._index.get(item_name)
        if vals is None:
            return None
        if source is not None:
            vals = [v for v in vals if v.source == source]
        if len(vals) == 1:
            return vals[0]
        # print("multiple with name", item_name)
        return None

    def check_uniqueness(self, prefer, index=None):
        if index is None:
            index = self._index
        for idx in index:
            if len(self._index[idx]) > 1:
                print(f"Multiple definitions found for {idx}:")
                for d in self._index[idx]:
                    print(f"  ->({d.source}){d.file}:{d.line_num}")


class Definitions:
    """
    Read a file and all of its includes.
    """

    def __init__(self, verbose=False):
        self.mol_defs = Index()
        self.atom_defs = Index()
        self.bond_types = Index()
        self.angle_types = Index()
        self.atoms = Index()
        self.bonds = Index()
        self.angles = Index()
        self.mols = Index()
        self.include_tree = defaultdict(list)
        self.include_lines = {}
        self.failed_includes = []
        self.loaded = set()
        self.verbose = verbose

    def get_molecules(self):
        """Get all instances of molecules. This essentially returns every occurence."""
        return [e for v in self.mols.values() for e in v]

    def get_molecule_names(self):
        return set(e.name for v in self.mols.values() for e in v)

    def get_atom_names(self):
        """Get all atom names"""
        return set(a.name for aas in self.atoms.values() for a in aas)

    def get_bond_names(self):
        """Get all atom names"""
        return set(a.name for aas in self.bonds.values() for a in aas)

    def get_angle_names(self):
        """Get all atom names"""
        return set(a.name for aas in self.angles.values() for a in aas)

    def load_itp_file(self, filename, source):
        """Load a file and all its dependencies. Use 'source' to track where original
        data came from, for example an actual topology file or manual addition"""
        if filename in self.loaded:
            return
        self.loaded.add(filename)
        self.load_includes(filename, source)
        fl = FileLoader(self, filename, source, self.verbose)
        fl.load_sections()

    def load_includes(self, filename, source):
        import os

        dirname = os.path.dirname(filename)
        with open(filename, "r") as f:
            lines = f.readlines()
            for idx, line in enumerate(lines):
                if line.startswith("#include"):
                    start, include_filename, end = line.split('"')
                    if dirname != "":
                        include_filename = dirname + "/" + include_filename
                    try:
                        self.include_lines[filename] = idx
                        self.include_tree[filename].append(include_filename)
                        self.load_itp_file(include_filename, source)
                    except Exception as e:
                        print(e)
                        self.failed_includes.append(include_filename)


class FileLoader:
    def __init__(self, defs, filename, source, verbose=False):
        self.defs = defs
        self.filename = filename
        self.line_num = 0
        self.mol_defs = 0
        self.atoms = 0
        self.bonds = 0
        self.mols = 0
        self.angles = 0
        self.atomtypes = 0
        self.bondtypes = 0
        self.angletypes = 0
        self.source = source
        self.verbose = verbose

        with open(filename, "r") as f:
            self.lines = f.readlines()

    def parse_moleculetype_atoms(self, mol_def):
        while True:
            line = self.next_line()
            if line is None:
                return
            atom_id = line.split()[0]
            atom_name = line.split()[1]
            charge = float(line.split()[6])
            self.defs.atoms.add(
                Atom(
                    atom_id,
                    atom_name,
                    charge,
                    self.filename,
                    self.line_num,
                    self.source,
                )
            )
            mol_def.atoms.add(
                Atom(
                    atom_id,
                    atom_name,
                    charge,
                    self.filename,
                    self.line_num,
                    self.source,
                )
            )
            self.atoms += 1

    def parse_moleculetype(self):
        molecule_name = self.next_line().split()[0]
        mol_def = self.defs.mol_defs.add(
            MolDef(molecule_name, self.filename, self.line_num, self.source)
        )
        self.mol_defs += 1
        while True:
            section = self.next_section()
            if section == "[atoms]":
                self.parse_moleculetype_atoms(mol_def)

            if section == "[bonds]":
                self.parse_moleculetype_bonds(mol_def)

            if section == "[angles]":
                self.parse_moleculetype_angles(mol_def)

            if section == "[moleculetype]":
                self.rewind()
                return None

            if section is None:
                return None

    def parse_moleculetype_bonds(self, mol_def):
        while True:
            line = self.next_line()
            if line is None:
                return
            bond = line.split()
            if len(bond) == 3:
                mol_def.bonds.add(
                    BondRef(bond[2], self.filename, self.line_num, self.source)
                )
                self.defs.bonds.add(
                    BondRef(bond[2], self.filename, self.line_num, self.source)
                )
                self.bonds += 1

    def parse_moleculetype_angles(self, mol_def):
        while True:
            line = self.next_line()
            if line is None:
                return
            bond = line.split()
            if len(bond) == 4:
                mol_def.angles.add(
                    AngleRef(bond[3], self.filename, self.line_num, self.source)
                )
                self.defs.angles.add(
                    AngleRef(bond[3], self.filename, self.line_num, self.source)
                )
                self.angles += 1

    def parse_atomtypes(self):
        while True:
            line = self.next_line()
            if line is None:
                return
            line = line.split()
            atom_name = line[0]
            self.defs.atom_defs.add(
                AtomDef(atom_name, self.filename, self.line_num, self.source)
            )
            self.atomtypes += 1

    def parse_bondtypes(self):
        while True:
            line = self.next_line()
            if line is None:
                return
            line = line.split()
            atom_name = line[0]
            self.defs.bond_types.add(
                BondType(atom_name, self.filename, self.line_num, self.source)
            )
            self.bondtypes += 1

    def parse_angletypes(self):
        while True:
            line = self.next_line()
            if line is None:
                return
            line = line.split()
            atom_name = line[0]
            self.defs.angle_types.add(
                AngleType(atom_name, self.filename, self.line_num, self.source)
            )
            self.angletypes += 1

    def parse_molecule(self):
        while True:
            line = self.next_line()
            if line is None:
                return
            line = line.split()
            mol_name = line[0]
            count = line[1]
            self.defs.mols.add(
                Mol(mol_name, int(count), self.filename, self.line_num, self.source)
            )
            self.mols += 1

    def next_line(self, check=True):
        while self.line_num < len(self.lines):
            line = self.lines[self.line_num]
            line = line.strip()
            line = line.split(";")[0]
            if check and len(line) > 0 and line[0] == "[":
                return None
            self.line_num += 1
            if line.startswith("#define"):
                return line.split("#define")[1]
            if len(line) == 0 or line[0] in ("#", ";"):
                continue
            else:
                return line

    def rewind(self):
        self.line_num -= 1

    def next_section(self):
        while self.line_num < len(self.lines):
            line = self.next_line(check=False)
            if line is None:
                return None
            line = line.strip()
            line = line.replace(" ", "")
            self.curr_line = line
            if len(line) > 0:
                if line[0] == "[":
                    return line
        return None

    def load_sections(self):
        if self.verbose:
            print("Loading", self.filename)
        curr_section = None
        with open(self.filename, "r") as f:
            self.lines = f.readlines()
            self.line_num = 0
        while self.line_num < len(self.lines):
            section = self.next_section()
            if section is None:
                return
            section = section.replace(" ", "")
            if section == "[moleculetype]":
                self.parse_moleculetype()
            if section == "[atomtypes]":
                self.parse_atomtypes()
            if section == "[molecules]":
                self.parse_molecule()
            if section == "[bondtypes]":
                self.parse_bondtypes()
            if section == "[angletypes]":
                self.parse_angletypes()
        if self.verbose:
            print("  molecule types:", self.mol_defs)
            print("  atoms :", self.atoms)
            print("  atomtypes :", self.atomtypes)
            print("  bondtypes :", self.bondtypes)
            print("  angletypes :", self.angletypes)
            print("  molecules :", self.mols)
