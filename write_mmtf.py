import mmtf
import mdtraj as mdt


def write_mmtf_file(traj: mdt.Trajectory) -> None:
    """ Writes an mdtraj trajectory to a mmtf file.
    """
    topology: mdt.Topology
    topology = traj.topology

    encoder = mmtf.MMTFEncoder()
    encoder.init_structure(
        total_num_bonds=topology.n_bonds,
        total_num_atoms=topology.n_atoms,
        total_num_chains=topology.n_chains,
        total_num_groups=topology.n_residues,
        total_num_models=1,
        structure_id="5ZMZ"
    )

    encoder.set_header_info(
        r_free=None,
        r_work=None,
        resolution=None,
        title=None,
        deposition_date=None,
        release_date=None,
        experimental_methods=None,
    )

    unit_cell = list(traj.unitcell_lengths[0, :]) + list(traj.unitcell_angles[0, :])
    assert len(unit_cell) == 6
    encoder.set_xtal_info(space_group="", unit_cell=unit_cell)

    # This sets the number of chains for a given model. Here we assume we have only
    # one model so, we only add one chain count
    encoder.set_model_info(model_id=None, chain_count=topology.n_chains)
    for chain in topology.chains:

        # encoder.set_entity_info(
        #     chain_indices=None,
        #     sequence=None,
        #     description=None,
        #     entity_type=None
        # )

        encoder.set_chain_info(
            chain_id=chain.index,
            chain_name="REPLACE ME",
            num_groups=chain.n_residues
        )

        for residue in chain.residues:

            encoder.set_group_info(
                group_name="REPLACE ME",
                group_number=
            )

    encoder.write_file("./temp.mmtf")


if __name__ == "__main__":
    trajectory = mdt.load("data/5zmz.pdb")
    write_mmtf_file(trajectory)
