import sys
import os
import glob
import itplz.mol_def
import itplz.tools.add_includes
import itplz.tools.charge
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
        topology = itplz.mol_def.Definitions(verbose=args.verbose)
        topology.load_itp_file(args.top, "top")
        mdls = load_external_mol_defs(
            args.directories, args.files, verbose=args.verbose
        )
        mdls.mols = topology.mols
        itplz.tools.charge.run(topology, mdls, args.verbose)

    if args.sub_name == "add-includes":
        topology = itplz.mol_def.Definitions(verbose=args.verbose)
        topology.load_itp_file(args.top, "top")
        mdls = load_external_mol_defs(
            args.directories, args.files, verbose=args.verbose
        )
        mdls.mols = topology.mols
        itplz.tools.add_includes.run(args.top, topology, mdls, args.out, args.verbose)
