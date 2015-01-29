=============
Hydrogen PDBs
=============

Upon reading in a PDB file, Rosetta protonates the PDB to match the Rosetta force field in a non-deterministic matter.
As this could affect the comparison of different protocols, the PDBs in this directory have been created by running them through Rosetta to protonate one time for all runs of a benchmark.
Rosetta is then run for then benchmark with the -no_optH flag set to true so that no alterations to hydrogens are made upon PDB loading, allowing a more direct comparison of different runs.
