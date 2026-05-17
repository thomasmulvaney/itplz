import sys


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def run(topology, mdls, verbose=False):
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
