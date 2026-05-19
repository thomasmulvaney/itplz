from collections import defaultdict


class Includes:
    def __init__(self):
        self.mapping = defaultdict(set)

    def add(self, file, reason):
        self.mapping[file].add(reason)


class Mapping:
    """Given two internal and external definitions, stores all definitions
    which are referenced, and all references which have no definitions in
    3 lists.

    - defined: All defininitions (internal or external) which are referenced.
    - external: All external definitions which are only external and referenced.
    - undefined: All definitions which are not defined.
    """

    def __init__(self, internal_defs, external_defs):
        self.internal_defs = internal_defs
        self.external_defs = external_defs
        self.defined = []
        self.undefined = []
        self.external = []

    def get_defs(self, refs):
        for r in refs.keys():
            if self.internal_defs.get(r):
                self.defined.append(self.internal_defs.get(r))
            elif self.external_defs.get(r):
                self.external.append(self.external_defs.get(r))
                self.defined.append(self.external_defs.get(r))
            else:
                self.undefined.append(r)


class DefMapping:
    """Stores the mappings for moltypes, atomtypes, bondtypes and angletypes."""

    def __init__(self, internal_defs, external_defs):
        self.mol_types = Mapping(internal_defs.mol_types, external_defs.mol_types)
        self.atom_types = Mapping(internal_defs.atom_types, external_defs.atom_types)
        self.bond_types = Mapping(internal_defs.bond_types, external_defs.bond_types)
        self.angle_types = Mapping(internal_defs.angle_types, external_defs.angle_types)

        # Find all defined molecules
        self.mol_types.get_defs(internal_defs.mols)

        # For every defined molecule, get all of the atoms, bonds and angle
        # types.
        for mol_type in self.mol_types.defined:
            self.atom_types.get_defs(mol_type.atoms)
            self.bond_types.get_defs(mol_type.bonds)
            self.angle_types.get_defs(mol_type.angles)

    def has_external(self):
        for t in [self.mol_types, self.atom_types, self.bond_types, self.angle_types]:
            if len(t.external) > 0:
                return True
        return False

    def has_undefined(self):
        for t in [self.mol_types, self.atom_types, self.bond_types, self.angle_types]:
            if len(t.undefined) > 0:
                return True
        return False

    def get_includes(self):
        """Returns includes"""
        incs = Includes()
        for atom_type in self.atom_types.defined:
            incs.add(atom_type.loc.file, "atomtypes")
        for bond_type in self.bond_types.defined:
            incs.add(bond_type.loc.file, "bondtypes")
        for angle_type in self.angle_types.defined:
            incs.add(angle_type.loc.file, "angletypes")
        for mol_type in self.mol_types.defined:
            incs.add(mol_type.loc.file, mol_type.name)

        return incs

    def get_ext_includes(self):
        """Returns any external includes"""
        incs = Includes()
        for atom_type in self.atom_types.external:
            incs.add(atom_type.loc.file, "atomtypes")
        for bond_type in self.bond_types.external:
            incs.add(bond_type.loc.file, "bondtypes")
        for angle_type in self.angle_types.external:
            incs.add(angle_type.loc.file, "angletypes")
        for mol_type in self.mol_types.external:
            incs.add(mol_type.loc.file, mol_type.name)

        return incs
