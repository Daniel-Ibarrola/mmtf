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

    unit_cell = [float(x) for x in traj.unitcell_lengths[0]]
    unit_cell += [float(x) for x in traj.unitcell_angles[0]]
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

        chain_name = chr(ord('A') + chain.index % 26)
        encoder.set_chain_info(
            chain_id=chain_name,
            chain_name="\x00",
            num_groups=chain.n_residues
        )

        for residue in chain.residues:

            encoder.set_group_info(
                group_name=residue.name,
                group_number=residue.index,
                insertion_code="\x00",
                group_type="",
                atom_count=residue.n_atoms,
                bond_count=0,  # This parameter is not used in mmtf-python
                single_letter_code=residue.code,
                sequence_index=-1,
                secondary_structure_type=-1,
            )

            for atom in residue.atoms:

                encoder.set_atom_info(
                    atom_name=atom.name,
                    serial_number=atom.serial,
                    alternative_location_id="\x00",
                    x=float(trajectory.xyz[0, atom.index, 0]),
                    y=float(trajectory.xyz[0, atom.index, 1]),
                    z=float(trajectory.xyz[0, atom.index, 2]),
                    occupancy=1.0,
                    temperature_factor=0,
                    element=atom.element.symbol,
                    charge=0,
                )

    for bond in topology.bonds:

        bond_order = bond.order
        if bond_order is None:
            bond_order = 1
        # Check if it's an intra-residue bond
        if bond.atom1.residue.index == bond.atom2.residue.index:

            encoder.current_group = encoder.group_list[bond.atom1.residue.index]

            residue_n_atoms = bond.atom1.residue.n_atoms
            atom_1_index = bond.atom1.index % residue_n_atoms
            atom_2_index = bond.atom2.index % residue_n_atoms
            # This method expects the atom index inside the residue or group
            encoder.set_group_bond(
                atom_index_one=atom_1_index,
                atom_index_two=atom_2_index,
                bond_order=bond_order
            )
        # If not it is an inter-residue bond
        else:
            # This method expects the atom index of the whole structure
            encoder.set_inter_group_bond(
                atom_index_one=bond.atom1.index,
                atom_index_two=bond.atom2.index,
                bond_order=bond_order
            )

    encoder.finalize_structure()
    encoder.write_file("./temp.mmtf")


if __name__ == "__main__":
    trajectory = mdt.load("data/5zmz.pdb")
    write_mmtf_file(trajectory)
