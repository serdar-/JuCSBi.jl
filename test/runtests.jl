using JuCSBi

atoms = parse_PDB("C:/Users/Serdar/PDB/1ggg.pdb.gz")
pdb = create_structure(atoms)
print(get_all_atoms(get_chain(pdb, "A"))[1])