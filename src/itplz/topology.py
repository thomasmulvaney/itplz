from itplz.mol_def import MolDefLoader

class Topology:
    def __init__(self, filename):
        self.filename = filename
        # All the molecule types read in the [ molecules ] field
        self.moleculetypes = set()
        # All the [ moleculetype ] defined in the file and includes.
        self.mol_defs = MolDefLoader()
        self.failed_includes = []
        self.mol_defs.load_itp_file(filename)

        with open(filename, 'r')  as f:
            lines = f.readlines()
            for  idx, line in enumerate(lines):
                with open(filename, 'r')  as f:
                    lines = f.readlines()
                    if line.startswith("[ molecules ]"):
                        self.read_molecules(lines[idx+1:])


    def read_molecules(self, lines):
        for line in lines:
            if line[0] == ';':
                continue
            if line[0] == '[':
                return
            self.moleculetypes.add(line.split()[0])
                
