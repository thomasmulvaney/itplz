import sys
import itplz.core


def check_type(type_map, name):
    t = len(type_map.defined)
    u = len(type_map.undefined)
    if len(type_map.undefined) == 0:
        if len(type_map.external) == 0:
            print(f"✔  {t}/{t} {name} types defined")
        else:
            print(
                f"✔  {t}/{t} {name} types defined, {len(type_map.external)} in external files"
            )
    else:
        print(f"✖  {t}/{t + u} {name} types defined")

    if len(type_map.undefined) > 0:
        print("   Undefined:")
        for m in type_map.undefined[0:5]:
            print(f"     {m:<20}")
        if len(type_map.undefined) > 5:
            print(f"     ... and {len(type_map.undefined) - 5} more")


def run(top, extra_mdls, args, verbose=False):
    def_map = itplz.core.DefMapping(top, extra_mdls)
    check_type(def_map.mol_types, "molecule")
    check_type(def_map.atom_types, "atom")
    check_type(def_map.bond_types, "bond")
    check_type(def_map.angle_types, "angle")

    if def_map.has_undefined():
        print()
        print("Some definitions are missing.")
        print(
            "Try adding additional files ITP using `--directories` or `--files` options"
        )
        sys.exit(2)

    if def_map.has_external():
        print()
        print("Some definitions are from externally included files:")
        for f in def_map.get_ext_includes().mapping.keys():
            print(f"  {f}")
        dir_arg = ""
        file_arg = ""
        if len(args.directories) > 0:
            dir_arg = f"--directories {' '.join(args.directories)}"

        if len(args.files) > 0:
            file_arg = f"--files {' '.join(args.files)}"
        print()
        print(
            f"Run `itplz add-includes --top {args.top} {dir_arg} {file_arg} --out path/to/new.top` to create a complete topology file."
        )
        sys.exit(1)

    print("Topology appears fine")
    sys.exit(0)
