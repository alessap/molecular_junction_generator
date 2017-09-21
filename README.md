# molecular_junction_generator

Generate a molecular junction (electrode-molecule-electrode) 
with FCC111 gold electrodes.

This is a simple script that takes a molecular structure XYZ file 
as input file and create a GEN (and XYZ) file of the molecule sandwiched 
between two gold FCC111 electrodes of size 4x4x6.

This script needs ASE installed in your system. It has been tested on macOSX.

The little guide below can be obtained running:
$ python generate_molecular_junction.py -h

Example 1:
$ python generate_molecular_junction.py molecule.xyz c r 0 1

using molecule 'molecule.xyz' as device region
it generates a cluster type molecular junction (c) with no pbc,
removes terminal hydrogen (r) on the end atoms 0 and 1.

Example 2:
$ python generate_molecular_junction.py molecule.xyz s nr

generate a supercell type molecular junction with pbc (s),
does not removes terminal hydrogen (nr) and automatically looks for
sulphur atoms as end atoms.
