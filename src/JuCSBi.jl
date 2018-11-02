module JuCSBi

using GZip
using Distances

PDB_DIR = joinpath(pwd(),"PDB") # Default PDB download directory
AW = Dict([("C",12.0), ("O",16.0), ("S",32.0), ("N", 14.0)]) # Atomic weights
AA_NAMES = Dict(["ALA"=>"A", "CYS"=>"C", "ASP"=>"D", "GLU"=>"E", "PHE"=>"F", # One letter aminoacid codings
"GLY"=>"G", "HIS"=>"H", "ILE"=>"I", "LYS"=>"K", "LEU"=>"L", "MET"=>"M",
"ASN"=>"N", "PRO"=>"P", "GLN"=>"Q", "ARG"=>"R", "SER"=>"S", "THR"=>"T",
"VAL"=>"V", "TRP"=>"W", "TYR"=>"Y"])

# Biological structure related types
"""
  Keeps the information for atoms provided by PDB file.
"""
struct Atom
  serial::Int32
  name::String
  residue_name::String
  chain::String
  residue_serial::Int32
  coords::Array{Float64,1}
  occupancy::Float64
  βfactor::Float64
  element::String
  charge::String
  alternative:: Bool
  Atom(serial, name, residue_name, chain, residue_serial, coords, occupancy, βfactor, element, charge, alternative) =
    new(serial, name, residue_name, chain, residue_serial, coords, occupancy, βfactor, element, charge, alternative)
end

"""
  Keeps the residue information provided by PDB file, contains an array of Atoms.
"""
struct Residue
  serial::Int32
  name::String
  chain::String
  atoms::Array{Atom,1}
  alternative::Bool
  Residue(serial, name, chain, atoms, alternative) = new(serial, name, chain, atoms, alternative)
end

"""
  Type that provides chain information, it holds an array of Atoms and Residues
  which constitute the chain.
"""
struct Chain
  name::String
  residues::Array{Residue,1}
  atoms::Array{Atom,1}
  Chain(name, residues, atoms) = new(name, residues, atoms)
end

"""
  βfactor that keeps model information which is parsed from PDB files with NMR
  structures. Currently not supported.
"""
struct Model
  serial::Int32
  chains::Array{Chain,1}
  Model(serial, chains) = new(serial, chains)
end

"""
  Keeps all structure related information that is provided by ATOM sections of
  PDB file. HETATM information is currently not supported.
"""
struct PDB
  name::String
  chains::Array{Chain,1}
  residues::Array{Residue,1}
  atoms::Array{Atom,1}
  PDB(name, chains, residues, atoms) = new(name, chains, residues, atoms)
end

# Methods for atom type filtering

isCA(atom::Atom) = atom.name == "CA"
isCB(atom::Atom) = atom.name == "CB"
isC(atom::Atom) = atom.name == "C"
isN(atom::Atom) = atom.name == "N"
isS(atom::Atom) = atom.name == "S"
isO(atom::Atom) = atom.name == "O"
ofN(atom::Atom) = contains(atom.name, "N")
ofC(atom::Atom) = contains(atom.name, "C")
ofO(atom::Atom) = contains(atom.name, "O")
ofS(atom::Atom) = contains(atom.name, "S")
on_backbone(atom::Atom) = isN(atom) | isC(atom) | isCA(atom)
on_sidechain(atom::Atom) = !isN(atom) | !isC(atom) | !isCA(atom)

# Methods for selecting structural elements (atoms, residues, chains)

get_chains(pdb::PDB, chain_id::String) = filter(c -> contains(chain_id, c.name), pdb.chains)
get_residue(pdb::PDB, chain_id::String, residue_serial::Int) = filter(res -> res.serial == residue_serial && res.chain == chain_id, pdb.residues)[1]
get_residue(chain::Chain, residue_serial::Int) = filter(res -> res.serial == residue_serial, chain.residues)[1]
get_residues_by_name(pdb::PDB, chain_id::String, residue_name::String) = filter(res -> res.name == residue_name, pdb.residues)
get_residues_by_name(chain::Chain, residue_name::String) = filter(res -> res.name == residue_name, chain.residues)
"Returns an array of atoms "
get_atoms(filter_func::Function, pdb::PDB) = filter(filter_func, pdb.atoms)
get_atoms(filter_func::Function, chain::Chain) = filter(filter_func, chain.atoms)

function get_atoms(filter_func::Function, chains::Array{Chain,1})
  atoms = Array(Atom,0)
  for chain in chains
    vcat(atoms, get_atoms(filter_func, chain))
  end
  atoms
end
get_atoms(filter_func::Function, residue::Residue) = filter(filter_func, residue.atoms)
get_atoms(filter_func::Function, atoms::Array{Atom,1}) = filter(filter_func, atoms)
"Returns an Array of alpha-carbon atoms."
get_CAs(pdb::PDB) = filter(isCA, pdb.atoms)
get_CAs(pdb::PDB, chain_id::String) = filter(atom -> atom.name == "CA" && atom.chain == chain_id, pdb.atoms)
get_CAs(chain::Chain) = filter(isCA, chain.atoms)
get_CAs(pdb::PDB, chain_ids::String) = filter(atom -> atom.name == "CA" && contains(chain_ids, atom.chain), pdb.atoms)
get_CAs(residue::Residue) = filter(isCA, residue.atoms)
get_CAs(residues::Array{Residue,1}) = map(get_CAs, residues)
get_CAs(atoms::Array{Atom,1}) = filter(atom -> atom.name =="CA", atoms)
get_coords(atoms::Array{Atom,1}) = map(atom -> atom.coords, atoms)
function get_coords(atoms::Array{Atom,1})
  onedim = map(atom -> atom.coords, atoms)
  coordinates = Array(Float64, length(onedim), 3)
  for i = 1:length(onedim)
    coordinates[i,:] = onedim[i]
  end
  coordinates
end

# PDB IO and parsing related functions

function download_PDB(PDBCode::String, downloadDir::String=PDB_DIR)
  try
    if !isdir(PDB_DIR)
      mkdir(PDB_DIR)
      println(string("Created PDB directory ", PDB_DIR))
    end
    filepath = joinpath(downloadDir,string(PDBCode, ".pdb.gz"))
    download(string("http://www.rcsb.org/pdb/files/" ,uppercase(PDBCode) ,".pdb.gz"), filepath)
    println(string(PDBCode, " is saved to ", downloadDir, "."))
  catch(e)
      println(string("There was a problem downloading ", PDBCode, "."))
      print(e)
  end
end

function read_PDB(PDBFile::String)
  pdbtext = nothing
  try
    if endswith(PDBFile, ".gz")
      stream = GZip.open(PDBFile, "r")
      pdbtext = readlines(stream)
    else
      stream = open(PDBFile, "r")
      pdbtext = readlines(stream)
    end
  catch(e)
    print(e)
  end
  pdbtext
end

function create_structure(atoms::Array{Atom,1}) # Add optional name argument to function
  residues = Residue[]
  chains = Chain[]
  residue_chain = unique([(atom.residue_serial,atom.chain) for atom in atoms])
  only_chains = unique([atom.chain for atom in atoms])
  for rc in residue_chain
    separate_res(atom) = atom.residue_serial == rc[1] && atom.chain == rc[2]
    residue_atoms = filter(separate_res, atoms)
    push!(residues, Residue(residue_atoms[1].residue_serial,residue_atoms[1].residue_name, residue_atoms[1].chain, residue_atoms, false))
  end
  for c in only_chains
    separate_chains(res) = res.chain == c
    separate_atom(atom) = atom.chain == c
    chain_residues = filter(separate_chains, residues)
    chain_atoms = filter(separate_atom, atoms)
    push!(chains, Chain(c, chain_residues, chain_atoms))
  end
  PDB("null", chains, residues, atoms) # Add optional name argument to function
end

function parse_PDB(PDBFile::String)
  atoms = Atom[]
  pdbtext = read_PDB(PDBFile)
  for line in pdbtext
    if startswith(line, "ATOM")
      serial = parse(Int,strip(line[7:11]))
      name = strip(line[13:16])
      if line[17] == " "
        alternative = false
      else
        alternative = true
      end
      residue_name = line[18:20]
      chain = string(line[22])
      residue_serial = parse(Int,strip(line[23:26]))
      coords = [parse(Float64,strip(line[31:38])), parse(Float64,strip(line[39:46])), parse(Float64,strip(line[47:54]))]
      occupancy = parse(Float64,strip(line[55:60]))
      βfactor = parse(Float64,strip(line[61:66]))
      element = strip(line[77:78])
      charge = strip(line[79:80])
      new_atom = Atom(serial, name, residue_name, chain, residue_serial, coords, occupancy, βfactor, element, charge, alternative)
      push!(atoms, new_atom)
    end
  end
  atoms
end

# Structure related methods

get_distance_vector(a1::Atom, a2::Atom) = a2.coords - a1.coords
get_distance(a1::Atom, a2::Atom) = norm(get_distance_vector(a1,a2))

function center_of_mass(atoms::Array{Atom,1})
  total_weight = 0.0
  total_weight_x_coords = 0.0
  for atom in atoms
    if ofC(atom)
      total_weight += AW["C"]
      total_weight_x_coords += AW["C"]*atom.coords
    elseif ofO(atom)
      total_weight += AW["O"]
      total_weight_x_coords += AW["O"]*atom.coords
    elseif ofN(atom)
      total_weight += AW["N"]
      total_weight_x_coords += AW["N"]*atom.coords
    elseif ofS(atom)
      total_weight += AW["S"]
      total_weight_x_coords += AW["S"]*atom.coords
    end
  end
  total_weight_x_coords/total_weight
end

get_sidechain_centroid(residue::Residue) = center_of_mass(get_atoms(on_sidechain, residue.atoms))

function get_sidechain_centroids(chain::Chain)
  centroids = zeros(length(chain.residues),3)
  for i = 1:len(chain.residues)
    centroids[i,:] = get_sidechain_centroids(chain.residues[i])
  end
  centroids
end

function get_backbone_length(chain::Chain)
  backbone_length = 0.0
  backbone_atoms = get_atoms(on_backbone, chain)
  for i = 1:length(backbone_atoms)-1
    backbone_length += get_distance(backbone_atoms[i], backbone_atoms[i+1])
  end
  backbone_length
end

function get_dihedral_angle(atom1::Atom, atom2::Atom, atom3::Atom, atom4::Atom)
  v1 = get_distance_vector(atom1, atom2)
  v2 = get_distance_vector(atom2, atom3)
  v3 = get_distance_vector(atom3, atom4)
  rad2deg(atan2(dot(cross(cross(v1,v2),cross(v2,v3)),b2/norm(b2)),dot(cross(v1,v2),cross(v2,v3))))
end

function get_dihedral_angle(atoms::Array{Atom,1})
  angle = NaN
  if length(atoms) == 4
    angle = get_dihedral_angle(atom[0], atom[1], atom[2], atom[3])
  else
    println("Atom array size should be 4!")
  end
  angle
end

function get_phi_angle(atomC1::Atom, atomN::Atom, atomCA::Atom, atomC2::Atom)
  if isC(atomC1) && isN(atomN) && isCA(atomCA) && isC(atomC2)
    get_dihedral_angle(atomC1,atomN,atomCA,atomC2)
  end
end

function get_psi_angle(atomN1::Atom, atomCA::Atom, atomC::Atom, atomN2::Atom)
  if isN(atomN1) && isCA(atomCA) && isC(atomC) && isN(atomN2)
    get_dihedral_angle(atomN1, atomCA, atomC, atomN2)
  end
end

function get_pairwise_distances(atoms1::Array{Atom,1}, atoms2::Array{Atom,1})
  coords1 = get_coords(atoms1)
  coords2 = get_coords(atoms2)
  pairwise(Euclidean(), transpose(coords1), transpose(coords2))
end

function get_interface_atoms(atoms1::Array{Atom,1}, atoms2::Array{Atom,1}, cutoff::Float64=5.0)
  n = length(atoms1)
  m = length(atoms2)
  pw_matrix = get_pairwise_distances(atoms1, atoms2)
  neighbors = ind2sub((n,m),find(pw_matrix .< cutoff))
  c1 = unique(neighbors[1])
  c2 = unique(neighbors[2])
  interface_atoms = Array(Atom,length(c1) + length(c2))
  for i = 1:length(c1)
        interface_atoms[i] = atoms1[c1[i]]
  end
  for j = 1:length(c2)
        interface_atoms[j+length(c1)] = atoms2[c2[j]]
  end
  interface_atoms
end

function get_interface_atoms(chain1::Chain, chain2::Chain, cutoff::Float64=5.0)
  atoms1 = chain1.atoms
  atoms2 = chain2.atoms
  get_interface_atoms(atoms1, atoms2, cutoff=cutoff)
end

function get_interface_atoms(chains1::Array{Chain,1}, chains2::Array{Chain,1}, cutoff::Float64=5.0)
  atoms1 = Array(Atom,0)
  atoms2 = Array(Atom,0)
  for chain in chains1
    vcat(atoms1, chain.atoms)
  end
  for chain in chains2
    vcat(atoms2, chain.atoms)
  end
  get_interface_atoms(atoms1, atoms2, cutoff=cutoff)
end

end 
