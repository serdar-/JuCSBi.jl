# JuCSBi.jl Computational Structural Biology Library for Julia

JuCSBi contains a set of tools for parsing PDB (Protein Data Bank) files and calculating structure related parameters. It includes functions which can calculate torsional and dihedral angles, centroids of residues and side chains. There are also functions where you can estimate positions of interface residues in case if it is a complex structure and calculate interatomic distances.

## Installation

In order to install the package, run the following code in Julia:

```
using Pkg
Pkg.clone("https://github.com/serdar-/JuCSBi.jl")
```

## Roadmap 

- Add normal mode analysis functions (ANM and GNM)