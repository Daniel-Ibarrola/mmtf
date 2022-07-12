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
    topology = traj.topology
    bonds_list = list(topology.bonds)
    print(bonds_list[0])


if __name__ == "__main__":
    # load_mmft_with_mdanalysis()
    load_gromacs_file()
