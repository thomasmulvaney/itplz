import sys
import os
import glob
import itplz.mol_def
import argparse
import pathlib
from collections import defaultdict


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


parser = argparse.ArgumentParser(
    prog="itplz",
    description="Pleasing tools for working with ITP and gromacs files",
    epilog="Good luck!",
)

parser.add_argument(
    "--verbose",
    action="store_true",
    required=False,
    help="Print extra information while running",
)

subparsers = parser.add_subparsers(dest="sub_name")
moltypes = subparsers.add_parser(
    "moltypes", help="print out all molecule types covered by the set of ITP files"
)
moltypes.add_argument(
    "--directories",
    nargs="*",
    default=[],
    help="Loads any file with .itp file in the directory and all subdirectories",
)
moltypes.add_argument("--files", nargs="*", default=[])


print_include = subparsers.add_parser(
    "print-includes", help="Generate an include string for a set of moltypes"
)
print_include.add_argument("--directories", nargs="*", default=[])
print_include.add_argument("--files", nargs="*", default=[])
print_include.add_argument("--moltypes", nargs="+", default=[])

add_include = subparsers.add_parser(
    "add-includes", help="Writes a topology file with additional includes"
)
add_include.add_argument(
    "--directories", nargs="*", default=[], help="Directories to find ITP files"
)
add_include.add_argument("--files", nargs="*", default=[], help="List of ITP files")
add_include.add_argument("--top", help="Input topology file")
add_include.add_argument("--out", help="Output topology file")


check_top = subparsers.add_parser(
    "check-top",
    help="Check a topology file for issues. Can be optionally run with extra ITP files to check outcome.",
)
check_top.add_argument("--directories", nargs="*", default=[])
check_top.add_argument("--files", nargs="*", default=[])
check_top.add_argument("--top")

charge = subparsers.add_parser(
    "charges",
    help="Check a topology file for issues. Can be optionally run with extra ITP files to check outcome.",
)
charge.add_argument("--directories", nargs="*", default=[])
charge.add_argument("--files", nargs="*", default=[])
charge.add_argument("--top")


def load_external_mol_defs2(mdls, directories, files):
    for d in directories:
        for filename in glob.glob(f"{d}/*.itp"):
            mdls.load_itp_file(filename, "ext")

    for f in files:
        mdls.load_itp_file(f, "ext")


def load_external_mol_defs(directories, files, verbose):
    mdls = itplz.mol_def.Definitions(verbose=verbose)
    for d in directories:
        for filename in glob.glob(f"{d}/*"):
            mdls.load_itp_file(filename, "ext")

    for f in files:
        mdls.load_itp_file(f, "ext")
    return mdls


def run_moltypes(directories=[], files=[]):
    mdls = load_external_mol_defs(directories, files)
    for ds in mdls.mol_defs.values():
        for d in ds:
            print(f"{d.name:<10}{d.file}:{d.line_num}")


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


def run_print_includes(directories=[], files=[], moltypes=[], verbose=False):
    mdls = load_external_mol_defs(directories, files, verbose=verbose)
    print("; Includes generated by itplz")
    create_includes(None, mdls, moltypes)


def run_charges(top, directories=[], files=[], verbose=False):
    topology = itplz.mol_def.Definitions(verbose=verbose)
    topology.load_itp_file(top, "top")
    mdls = load_external_mol_defs(directories, files, verbose)
    mdls.mols = topology.mols
    total = 0
    eprint(f"{'Molecule':12} {'Charge':>10} {'Count':>10} {'Tot. charge':>13}")
    eprint("-" * 48)
    undefined = False
    for molecule in topology.get_molecules():
        mol_def = topology.mol_defs.get(molecule.name)
        if mol_def is None:
            mol_def = mdls.mol_defs.get(molecule.name)
        if mol_def is None:
            undefined = True
            eprint(f"{molecule.name:12} {'-':>10} {'-':>10} {'-':>13} ; undefined!")
        else:
            q = mol_def.charge()
            cnt = molecule.count
            eprint(
                f"{mol_def.name:12} {q:>10} {cnt:>10} {cnt * q:>13} ; ({mol_def.source}){mol_def.file}:{mol_def.line_num}"
            )
            total += q * cnt
    if undefined:
        eprint("Total Charge could not be determined")
        sys.exit(1)
    eprint("Total Charge:", total)
    print(total)
    sys.exit(0)


def run_add_includes(top, out_file, directories=[], files=[], verbose=False):
    topology = itplz.mol_def.Definitions(verbose=verbose)
    topology.load_itp_file(top, "top")
    mdls = load_external_mol_defs(directories, files, verbose=verbose)
    mdls.mols = topology.mols
    try:
        last_include = max(topology.include_lines.values())
    except Exception as e:
        last_include = 0
    incs = create_includes(topology, mdls, topology.get_molecule_names(), out_file)

    with open(out_file, "w") as f:
        for inc in incs.mapping:
            i_path = os.path.relpath(inc, os.path.dirname(out_file))
            f.write(f'#include "{i_path}" ; {", ".join(incs.mapping[inc])}\n')

        with open(top, "r") as top_f:
            for line in top_f.readlines():
                if not line.startswith("#include"):
                    f.write(line)


def run_check_top(top_file, extra_directories=[], extra_files=[], verbose=False):
    top = itplz.mol_def.Definitions(verbose=verbose)
    extra_mdls = itplz.mol_def.Definitions(verbose=verbose)
    top.load_itp_file(top_file, "top")
    load_external_mol_defs2(extra_mdls, extra_directories, extra_files)
    not_in_top = []
    not_in_extras = []
    print("✔  Molecule types which have moleculartype definitions already")
    for mt in top.get_molecule_names():
        if top.mol_defs.get(mt, source="top"):
            include = top.mol_defs.get(mt, source="top")
            print(f"{mt:<10} ({include.source}){include.file}:{include.line_num}")
        else:
            not_in_top.append(mt)
    print("# Bonds with out Bond types")
    for b in extra_mdls.get_bond_names():
        print(b, extra_mdls.bond_types.get(b))
    print()
    print("# Molecules which have types defined in extra itp files")
    for mt in not_in_top:
        if extra_mdls.mol_defs.get(mt) is not None:
            include = extra_mdls.mol_defs.get(mt)
            print(f"{mt:<10} ({include.source}){include.file}:{include.line_num}")
        else:
            not_in_extras.append(mt)
    print()
    print("# Molecules which lack a definition")
    for mt in not_in_extras:
        print(f"{mt}")

    print()
    print("# Includes which could not be loaded")
    for inc in top.failed_includes:
        print(f"{inc}")

    print()
    print("# Checking for multiple molecule definitions")
    extra_mdls.mol_defs.check_uniqueness(extra_mdls.mols.keys())

    print()
    print("# Checking for multiple atom definitions")
    extra_mdls.atom_defs.check_uniqueness(extra_mdls.atoms.keys())

    print()
    print("# Atoms without definitions")
    defined_atoms = set()
    # Map the forcefields to atoms which they define
    extra_defined_atoms = set()
    undefined_atoms = set()
    for atom in top.get_atom_names():
        if atom in top.atom_defs.keys():
            defined_atoms.add(atom)
        elif atom in extra_mdls.atom_defs.keys():
            # Really there could be multiple FF matches.
            ff = extra_mdls.atom_defs.get(atom)
            if ff is not None:
                extra_defined_atoms.add(ff)
        else:
            undefined_atoms.add(atom)
    print(
        f"{len(defined_atoms)} have been defined, you are probably missing the forcefield ITP."
    )
    includes = set()
    for a in defined_atoms:
        w = extra_mdls.atom_defs.get(a.name)
        includes.add((w.source, w.file))
    for source, file in includes:
        print(f"  ({source}){file}")
    print(f"Defined in external force field:")
    for ff in extra_defined_atoms:
        print(f"  {ff.name} ({ff.source}){ff.file}")
    print(
        f"{len(undefined_atoms)} have not been defined, you are probably missing the forcefield ITP."
    )
    for a in undefined_atoms:
        print(a)


def main():
    args = parser.parse_args()
    if args.sub_name == "moltypes":
        run_moltypes(args.directories, args.files)

    if args.sub_name == "print-includes":
        run_print_includes(args.directories, args.files, args.moltypes)

    if args.sub_name == "check-top":
        run_check_top(args.top, args.directories, args.files, args.verbose)

    if args.sub_name == "charges":
        run_charges(args.top, args.directories, args.files, args.verbose)

    if args.sub_name == "add-includes":
        run_add_includes(args.top, args.out, args.directories, args.files, args.verbose)
