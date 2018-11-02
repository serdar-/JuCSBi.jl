using JuCSBi

atoms = JuCSBi.parse_PDB("C:/Users/Serdar/Documents/GitHub/JuCSBi.jl/test/PDB/1ggg.pdb.gz")
pdb = JuCSBi.create_structure(atoms)
print(JuCSBi.get_all_atoms(JuCSBi.get_chain(pdb, "A"))[1])