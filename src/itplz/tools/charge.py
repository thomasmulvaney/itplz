import os
import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run(topology, mdls, verbose=False):
    total = 0
    eprint(f"{'Molecule':12} {'Charge':>10} {'Count':>10} {'Tot. charge':>13}")
    eprint("-" * 48)
    undefined = False
    for molecule in topology.get_molecules():
        mol_type = topology.mol_types.get(molecule.name)
        if mol_type is None:
            mol_type = mdls.mol_types.get(molecule.name)
        if mol_type is None:
            undefined = True
            eprint(f"{molecule.name:12} {'-':>10} {'-':>10} {'-':>13} ; undefined!")
        else:
            q = mol_type.charge()
            cnt = molecule.count
            eprint(
                f"{mol_type.name:12} {q:>10} {cnt:>10} {cnt * q:>13} ; {os.path.normpath(mol_type.loc.file)}:{mol_type.loc.line_num}"
            )
            total += q * cnt
    if undefined:
        eprint("Total Charge could not be determined")
        sys.exit(1)
    eprint("Total Charge:", total)
    print(total)
    sys.exit(0)
