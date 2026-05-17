import os
import glob
import itplz.core
import itplz.mol_def
import itplz.tools.add_includes
import itplz.tools.charge
import itplz.tools.check
import argparse

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

subparsers = parser.add_subparsers(title="subcommand", dest="sub_name")

add_include = subparsers.add_parser(
    "add-includes", help="Writes a topology file with additional includes"
)
add_include.add_argument(
    "--directories", nargs="*", default=[], help="Directories to find ITP files"
)
add_include.add_argument("--files", nargs="*", default=[], help="List of ITP files")
add_include.add_argument("--top", help="Input topology file")
add_include.add_argument("--out", required=True, type=str, help="Output topology file")


check_top = subparsers.add_parser(
    "check",
    help="Check a topology file to see if all types are defined.  Extra ITP files can be supplied with --directories and --files.",
)
check_top.add_argument("--directories", nargs="*", default=[])
check_top.add_argument("--files", nargs="*", default=[])
check_top.add_argument("--top")

charge = subparsers.add_parser(
    "charge",
    help="Check a topology file for issues. Can be optionally run with extra ITP files to check outcome.",
)
charge.add_argument("--directories", nargs="*", default=[])
charge.add_argument("--files", nargs="*", default=[])
charge.add_argument("--top")


def load_external(args):
    mdls = itplz.mol_def.Definitions(verbose=args.verbose)
    for d in args.directories:
        for filename in glob.glob(f"{d}/*.itp"):
            mdls.load_itp_file(filename, "ext")

    for f in args.files:
        mdls.load_itp_file(f, "ext")
    return mdls


def load_topology(args):
    topology = itplz.mol_def.Definitions(verbose=args.verbose)
    topology.load_itp_file(args.top, "top")
    return topology


def main():
    args = parser.parse_args()
    if args.sub_name is None:
        parser.print_help()

    if args.sub_name == "check":
        topology = load_topology(args)
        mdls = load_external(args)
        itplz.tools.check.run(topology, mdls, args, args.verbose)

    if args.sub_name == "charge":
        topology = load_topology(args)
        mdls = load_external(args)
        itplz.tools.charge.run(topology, mdls, args.verbose)

    if args.sub_name == "add-includes":
        topology = load_topology(args)
        mdls = load_external(args)
        itplz.tools.add_includes.run(args.top, topology, mdls, args.out, args.verbose)
