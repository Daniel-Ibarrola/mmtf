import mdtraj as mdt
import mdtraj.core.element as mdt_element

import copy
import os
import pprint
import xml.etree.ElementTree as ETree


def create_topology() -> mdt.Topology:
    topo = mdt.Topology()

    # Adding a new chain is pretty straightforward
    chain = topo.add_chain()
    # Adding a residue
    residue = topo.add_residue(name="ASP", chain=chain, resSeq=1)
    # Adding an atom
    element = mdt_element.get_by_symbol("C")
    topo.add_atom(name="CA",
                  element=element,
                  residue=residue)

    element = mdt_element.get_by_symbol("N")
    topo.add_atom(name="NZ",
                  element=element,
                  residue=residue)

    # Adding a bond
    topo.add_bond(atom1=topo.atom(0),
                  atom2=topo.atom(1)
                  )
    return topo


def load_name_replacement_tables():
    """Load the list of atom and residue name replacements."""

    residue_name_replacements = {}
    atom_name_replacements = {}

    tree = ETree.parse(os.path.join(os.path.dirname(__file__), 'data', 'pdbNames.xml'))
    all_residues = {}
    protein_residues = {}
    nucleic_acid_residues = {}
    for residue in tree.getroot().findall('Residue'):
        name = residue.attrib['name']
        if name == 'All':
            parse_residue_atoms(residue, all_residues)
        elif name == 'Protein':
            parse_residue_atoms(residue, protein_residues)
        elif name == 'Nucleic':
            parse_residue_atoms(residue, nucleic_acid_residues)
    for atom in all_residues:
        protein_residues[atom] = all_residues[atom]
        nucleic_acid_residues[atom] = all_residues[atom]
    for residue in tree.getroot().findall('Residue'):
        name = residue.attrib['name']
        for id_ in residue.attrib:
            if id_ == 'name' or id_.startswith('alt'):
                residue_name_replacements[residue.attrib[id_]] = name
        if 'type' not in residue.attrib:
            atoms = copy.copy(all_residues)
        elif residue.attrib['type'] == 'Protein':
            atoms = copy.copy(protein_residues)
        elif residue.attrib['type'] == 'Nucleic':
            atoms = copy.copy(nucleic_acid_residues)
        else:
            atoms = copy.copy(all_residues)
        parse_residue_atoms(residue, atoms)
        atom_name_replacements[name] = atoms

    return residue_name_replacements, atom_name_replacements


def parse_residue_atoms(residue, dictionary):
    for atom in residue.findall('Atom'):
        name = atom.attrib['name']
        for id_ in atom.attrib:
            dictionary[atom.attrib[id_]] = name


def print_pdb_residue_name_replacements():

    residue_replacements, atom_replacements = load_name_replacement_tables()
    print("Residue replacements")
    print(f"Number of residue replacements {len(residue_replacements)}")
    pprint.pprint(residue_replacements)
    print("\nAtom replacements")
    print(f"Atom replacements length {len(atom_replacements)}")
    pprint.pprint(atom_replacements)


if __name__ == "__main__":

    topology = create_topology()

    print(f"Num chains:  {topology.n_chains}")
    print(f"Num residues:  {topology.n_residues}")
    print(f"Num atoms:  {topology.n_atoms}")
    print(f"Num bonds:  {topology.n_bonds}")

    print_pdb_residue_name_replacements()
