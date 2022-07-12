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

    atom_index = 0
    groups_per_chain = []

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

            bond_count = []

            encoder.set_group_info(
                group_name=residue.name,
                group_number=residue.index,
                insertion_code="\x00",
                group_type="",
                atom_count=residue.n_atoms,
                bond_count=bond_count,
                single_letter_code=residue.code,
                sequence_index=None,
                secondary_structure_type=-1,
            )

            for atom in residue.atoms:

                encoder.set_atom_info(
                    atom_name=atom.name,
                    serial_number=atom.serial,
                    alternative_location_id="\x00",
                    x=trajectory.xyz[0, atom_index, 1],
                    y=trajectory.xyz[0, atom_index, 2],
                    z=trajectory.xyz[0, atom_index, 3],
                    occupancy=None,
                    temperature_factor=None,
                    element=atom.element.name,
                    charge=0,
                )

    encoder.groups_per_chain = groups_per_chain

    encoder.finalize_structure()
    encoder.write_file("./temp.mmtf")


if __name__ == "__main__":
    trajectory = mdt.load("data/5zmz.pdb")
    write_mmtf_file(trajectory)
