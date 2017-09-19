#! /usr/bin/env python

import os
from sys import argv
import numpy as np
from ase import Atoms
from ase.io import read,write
from ase.lattice.surface import fcc111
from gpaw.transport.tools import sort_atoms

def find_nearest(atom, atoms):
    #create empty list
    dists = []
    #define coordinates for atom
    coord = atom.position
    #loop over atoms and append distance to list
    for a in atoms:
        if a.symbol=='H':
            v = a.position -coord
            dist = np.sqrt(v[0]**2+v[1]**2+v[2]**2)
            dists.append(dist)
            continue
        dists.append(1.7)
    #make list to np.array
    dists= np.array(dists)
    #find index with shortest distance using argmin
    index = np.argmin(dists)
    if dists[index]>1.7: # atom closer than C
        print 'found no H bonded to S, deleted something anyway.consider revising'
    return index

#mol=read('mol.xyz')
#natoms=len(mol)

if len(argv) > 3:
    print("reading end atoms from lunch")
    endatom1=int(argv[2])
    endatom2=int(argv[3])

#create gold electrode
a=4
pbc=bool(True)

b=a
c=6
n=0
c=c+n
lelectrode = fcc111('Au', size=(a,b,c), vacuum=5.0,orthogonal=True) 
# Au-S 2.72AA / Au-Au 2.88AA
#deltag           2.35558909829 # distance between two adjacent gold layers
#delta s_g        2.14441090171 # perpendicular distance between endatom and electrode surface

# read the molecule file
molecule = read('mol.xyz')
#place electrode close to the origin of axis
n=(lelectrode[0].position)
n = - n                                 # put atom 1 at the origin of axis
lelectrode.translate(n)
natoms = len(molecule)

try: 
    endatom1 
except NameError:
    endatom1 = -1
    for atom in molecule:
        if atom.symbol == 'S': # is your end atoms are sulphur atoms and there are only 2 S atoms in the molecule
            if endatom1 < 0:
                endatom1=atom.index
            else:
                endatom2=atom.index

#remove terminal hydrogens
#list =[]
#for atom in molecule:
#    if atom.symbol =='S':
#        list.append(atom.index)
#del molecule[find_nearest(molecule[list[1]],molecule)]
#del molecule[find_nearest(molecule[list[0]],molecule)]
#natoms = len(molecule)

# align the sulphur atoms along the z axis
po=(molecule[endatom1].position)
lo=(molecule[endatom2].position)
v=lo-po
z=[0,0,1]
molecule.rotate(v,z)
molecule.rotate('z',np.pi+np.pi/4.0)

# place the molecule on top of a hollow fcc site for axa face
if a == 6:
    q=( lo + [0,0,4.5] - lelectrode[(b*a*2)+(6*4)].position )
if a == 3:
    q=( lo + [0,0,4.5] - lelectrode[(a*b*2)+(4)].position )
if a == 4:
    q=( lo + [0,0,4.5] - lelectrode[(b*a*2)+(9)].position )

q = - q
molecule.translate(q)
#define distance betweet gold layers
deltag=abs(lelectrode[0].position[2] - lelectrode[(a*b)+1].position[2])
#define end to end distance
v=molecule[endatom2].position-molecule[endatom1].position
lenght_molecule=np.sqrt(v[0]**2+v[1]**2+v[2]**2)
molecule.translate([0,0,-(deltag)])
molecule.extend(lelectrode)
# make a list of indexes relative to gold atoms
au_list =[]
for atom in molecule:
    if atom.symbol =='Au':
      au_list.append(atom.index)

# number of gold atoms on a single layer
ng=int(a*b)
# distance along the z axis between one layer and the adjacent one:
deltag=abs(lelectrode[0].position[2] - lelectrode[(a*b)+1].position[2])
v=molecule[endatom2].position-molecule[natoms+2].position
deltas_g=abs(v[2])

for x in range(c):
    ly=[]
    ly=au_list[(x*a*b):(a*b*(x+1))]
    layer=molecule[ly]
    layer.translate([0,0,-(abs(lenght_molecule)+abs(deltas_g*2)+abs(deltag*(2*x)))])
    molecule.extend(layer)

if argv[1] == "s" or argv[1] == "S":
    molecule.set_pbc(True)
if argv[1] == "c" or argv[1] == "C":
    molecule.set_pbc(False)
n=(molecule[natoms+a*b*c-a*b].position)
n=-n
molecule.translate(n)

if argv[1] == "s" or argv[1] == "S":
    lele = fcc111('Au', size=(a,b,18), vacuum=deltas_g+lenght_molecule,orthogonal=True)
    cell=lele.get_cell()
    molecule.set_cell([cell[0],cell[1],[  0. , 0. , molecule[len(molecule)-1].position[2]-20.0]])
    molecule.set_pbc((True, True, True))
    write('moljunct_S.gen',molecule)
    write('moljunct_S.xyz',molecule)
    write('xyz_moljunct_S.xyz',molecule)
write('moljunct.gen',molecule)
write('xyz_moljunct.xyz',molecule)


# printing device, source and drain atom indeces
print 1, natoms 
print natoms+1 , natoms+a*b*6
print natoms+a*b*6+1, len(molecule)

exit()
