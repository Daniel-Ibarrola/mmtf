import to_mdtraj as conv
import MDAnalysis as mda
import mdtraj as mdt
import mmtf


def load_mmft_with_custom_module():
    traj = conv.load_mmtf("1sux.mmtf")
    print(traj)


def load_mmft_with_mdanalysis():
    universe = mda.Universe("1sux.mmtf")
    print(universe)
    print(f"Number of residues: {len(universe.residues)}")
    print(universe.atoms.positions.shape)
    print(f"Number of models: {len(universe.models)}")


def load_gromacs_file():
    traj = mdt.load("md_1u19.gro")
    assert isinstance(traj, mdt.Trajectory)
    print(traj)
    print(traj.n_atoms)


def parse_with_mmtf():
    decoder = mmtf.parse("1sux.mmtf")
    # With decoder.group_type_list we can obtain a list of the group
    # to which each resiude belongs to. Groups contain information
    # about atom names, bonds, etc

    # Group type list is a list of dictionaries where each entry
    # has the following structure
    # {'groupName': 'ALA',
    #  'atomNameList': ['N', 'CA', 'C', 'O', 'CB'],
    #  'elementList': ['N', 'C', 'C', 'O', 'C'],
    #  'bondOrderList': [1, 1, 2, 1],
    #  'bondAtomList': [1, 0, 2, 1, 3, 2, 4, 1],
    #  'formalChargeList': [0, 0, 0, 0, 0],
    #  'singleLetterCode': 'A',
    #  'chemCompType': 'L-PEPTIDE LINKING'}


if __name__ == "__main__":
    load_mmft_with_mdanalysis()
    # load_gromacs_file()
