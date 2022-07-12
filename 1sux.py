import mmtf
import mdtraj as mdt


def load_1sux_with_mmtf(decoder):

    print(decoder.structure_id)
    print(f"Num models:  {decoder.num_models}")
    print(f"Num chains:  {decoder.num_chains}")
    print(f"Num groups:  {decoder.num_groups}")
    print(f"Num atoms:  {decoder.num_atoms}")
    print(f"Num bonds:  {decoder.num_bonds}")


def traverse_mmtf(decoder):

    model_index = 0
    atom_index = 0
    chain_index = 0
    group_index = 0

    for model_chain_count in decoder.chains_per_model:
        print(f"Model: {model_index}")
        for chain in range(model_chain_count):
            # Note: chain names in chain_name_list are not unique
            chain_name = decoder.chain_name_list[chain_index]
            # Number of groups per chain
            chain_group_count = decoder.groups_per_chain[chain_index]
            print(f"Chain {chain_name} with index {chain} has {chain_group_count} group(s)")
            for ii in range(chain_group_count):
                group = decoder.group_list[decoder.group_type_list[group_index]]
                group_name = group["groupName"]
                # Note: ions and water have '?' for singel letter code
                group_single_letter_code = group["singleLetterCode"]
                print(f"Group name is {group_name} ({group_single_letter_code})")

                atom_offset = atom_index
                group_bond_count = int(len(group["bondAtomList"]) / 2)

                # Traverse bonds between atoms inside groups
                # print("Bonds")
                # for ii in range(group_bond_count):
                #     print(atom_offset + group["bondAtomList"][ii * 2])
                #     print(atom_offset + group["bondAtomList"][ii * 2 + 1])

                group_atom_count = len(group["atomNameList"])
                for jj in range(group_atom_count):
                    coords = [
                        decoder.x_coord_list[atom_index],
                        decoder.y_coord_list[atom_index],
                        decoder.z_coord_list[atom_index],
                              ]
                    atom_name = group["atomNameList"][jj]
                    element = group["elementList"][jj]
                    print(f"Atom {atom_name}. Element {element}")
                    print(coords)
                    atom_index += 1
                group_index += 1
            chain_index += 1
        model_index += 1


def load_1sux_with_mdtraj():
    traj = mdt.load("1sux.pdb")
    topology = traj.topology
    assert isinstance(topology, mdt.Topology)

    print(f"Num chains:  {topology.n_chains}")
    print(f"Num residues:  {topology.n_residues}")
    print(f"Num atoms:  {topology.n_atoms}")
    print(f"Num bonds:  {topology.n_bonds}")

    print("A chain:", topology.chain(0), type(topology.chain(0)))
    print("A residue:", topology.residue(0), type(topology.residue(0)))
    print("An atom:", topology.atom(0), type(topology.atom(0)))


if __name__ == "__main__":

    print("Loading mmtf")
    mmtf_decoder = mmtf.parse("1sux.mmtf")
    load_1sux_with_mmtf(mmtf_decoder)
    traverse_mmtf(mmtf_decoder)

    # print("\nLoading mdtraj")
    # load_1sux_with_mdtraj()
