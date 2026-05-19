from collections import defaultdict


class Loc:
    """A location of a definition"""

    __slots__ = ("file", "line_num")

    def __init__(self, file, line_num):
        self.file = file
        self.line_num = line_num


class Definition:
    """
    Instance of a definition of something.
    """

    __slots__ = ("name", "loc")

    def __init__(self, name: str, loc: Loc):
        self.name = name
        self.loc = loc

    def key(self):
        return self.name


class MolType(Definition):
    def __init__(self, name: str, loc: Loc):
        self.name = name
        self.loc = loc
        self.atoms = RefIndex()
        self.bonds = RefIndex()
        self.angles = RefIndex()

    def get_atom_names(self):
        return [a.name for aas in self.atoms.values() for a in aas]

    def get_bond_names(self):
        return [a.name for aas in self.bonds.values() for a in aas]

    def get_angle_names(self):
        return [a.name for aas in self.angles.values() for a in aas]

    def charge(self):
        total = 0
        for atoms in self.atoms.values():
            for a in atoms:
                total += a.charge
        return total


class AtomType(Definition):
    pass


class BondType(Definition):
    pass


class AngleType(Definition):
    pass


class BondRef(Definition):
    """A bond that depends on a BondType definition.

    Looks like:

       int int ref

    We just store the ref for now.
    """

    pass


class AngleRef(Definition):
    """A bond that depends on a BondType definition.

    Looks like:

       int int int ref

    We just store the ref for now.
    """

    pass


class Mol(Definition):
    """A molecule, whose name references a MolType."""

    __slots__ = ("molecule", "count", "loc")

    def __init__(self, name: str, count: int, loc: Loc):
        self.name = name
        self.count = count
        self.loc = loc


class Atom(Definition):
    __slots__ = ("id", "name", "charge", "loc")

    def __init__(self, id: int, name: str, charge: float, loc: Loc):
        self.id = id
        self.name = name
        self.charge = charge
        self.loc = loc

    def key(self):
        return self.name


class Index:
    def get(self, item):
        return self._index.get(item)

    def keys(self):
        return self._index.keys()

    def values(self):
        return self._index.values()


class DefIndex(Index):
    """Store definitions. Definition may only occur once."""

    def __init__(self):
        self._index = {}

    def add(self, item):
        if self._index.get(item.key()):
            f = self._index.get(item.key())
            print(
                "already defined!",
                item.name,
                item.file,
                item.line_num,
                f.name,
                f.file,
                f.line_num,
            )
        else:
            self._index[item.key()] = item
        return item


class RefIndex(Index):
    """Store references to definitions. References may occur multiple times.

    References do not need to have definitions. We check this elsewhere.
    """

    def __init__(self):
        self._index = defaultdict(list)

    def add(self, item):
        self._index[item.key()].append(item)


class Definitions:
    """
    Read a file and all of its includes.
    """

    def __init__(self, verbose=False):
        self.mol_types = DefIndex()
        self.atom_types = DefIndex()
        self.bond_types = DefIndex()
        self.angle_types = DefIndex()
        self.mols = RefIndex()
        self.include_tree = defaultdict(list)
        self.include_lines = {}
        self.failed_includes = []
        self.loaded = set()
        self.verbose = verbose

    def add_mol_type(self, obj):
        return self.mol_types.add(obj)

    def add_atom_type(self, obj):
        self.atom_types.add(obj)

    def add_bond_type(self, obj):
        self.bond_types.add(obj)

    def add_angle_type(self, obj):
        self.angle_types.add(obj)

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

    def load_itp_file(self, filename):
        """Load a file and all its dependencies."""
        if filename in self.loaded:
            return
        self.loaded.add(filename)
        self.load_includes(filename)
        fl = FileLoader(self, filename, self.verbose)
        fl.load_sections()

    def load_includes(self, filename):
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
                        self.load_itp_file(include_filename)
                    except Exception:
                        self.failed_includes.append(include_filename)


class FileLoader:
    def __init__(self, defs, filename, verbose=False):
        self.defs = defs
        self.filename = filename
        self.line_num = 0
        self.mol_types = 0
        self.atoms = 0
        self.bonds = 0
        self.mols = 0
        self.angles = 0
        self.atomtypes = 0
        self.bondtypes = 0
        self.angletypes = 0
        self.verbose = verbose

        with open(filename, "r") as f:
            self.lines = f.readlines()

    def loc(self):
        return Loc(self.filename, self.line_num)

    def parse_moleculetype_atoms(self, mol_type):
        while True:
            line = self.next_line()
            if line is None:
                return
            atom_id = line.split()[0]
            atom_name = line.split()[1]
            charge = float(line.split()[6])
            mol_type.atoms.add(
                Atom(
                    atom_id,
                    atom_name,
                    charge,
                    self.loc(),
                )
            )
            self.atoms += 1

    def parse_moleculetype(self):
        molecule_name = self.next_line().split()[0]
        mol_type = self.defs.add_mol_type(MolType(molecule_name, self.loc()))
        self.mol_types += 1
        while True:
            section = self.next_section()
            if section == "[atoms]":
                self.parse_moleculetype_atoms(mol_type)

            if section == "[bonds]":
                self.parse_moleculetype_bonds(mol_type)

            if section == "[angles]":
                self.parse_moleculetype_angles(mol_type)

            if section == "[moleculetype]":
                self.rewind()
                return None

            if section is None:
                return None

    def iter_lines(self):
        while True:
            line = self.next_line()
            if line is not None:
                yield (line, self.loc())
            else:
                return

    def parse_moleculetype_bonds(self, mol_type):
        for line, loc in self.iter_lines():
            bond = line.split()
            if len(bond) == 3:
                mol_type.bonds.add(BondRef(bond[2], loc))
                self.bonds += 1

    def parse_moleculetype_angles(self, mol_type):
        for line, loc in self.iter_lines():
            bond = line.split()
            if len(bond) == 4:
                mol_type.angles.add(AngleRef(bond[3], loc))

    def parse_atomtypes(self):
        for line, loc in self.iter_lines():
            line = line.split()
            atom_name = line[0]
            self.defs.atom_types.add(AtomType(atom_name, loc))
            self.atomtypes += 1

    def parse_bondtypes(self):
        for line, loc in self.iter_lines():
            line = line.split()
            atom_name = line[0]
            self.defs.bond_types.add(BondType(atom_name, loc))
            self.bondtypes += 1

    def parse_angletypes(self):
        for line, loc in self.iter_lines():
            line = line.split()
            atom_name = line[0]
            self.defs.angle_types.add(AngleType(atom_name, loc))
            self.angletypes += 1

    def parse_molecule(self):
        for line, loc in self.iter_lines():
            line = line.split()
            mol_name = line[0]
            count = line[1]
            self.defs.mols.add(Mol(mol_name, int(count), self.loc()))
            self.mols += 1

    def next_line(self, check=True):
        """Reads lines in a section. At the end of a section returns None"""
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
            print("  molecule types:", self.mol_types)
            print("  atoms :", self.atoms)
            print("  atomtypes :", self.atomtypes)
            print("  bondtypes :", self.bondtypes)
            print("  angletypes :", self.angletypes)
            print("  molecules :", self.mols)
