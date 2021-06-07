import CDPL.ForceField as ForceField
import numpy as np

from src.utils.Protein import Protein

# load protein
p = Protein()
p.fromFile('../data/../data/benchmark5.5/structures/1A2K_l_b.pdb')
p.prepare()  # calculates H-Atoms!
coords = p.getCoordinates().toArray(False)

randomSurfacePoints = np.random.randint(0, 100, size=(100, 3))

# prepare protein for force field calculations
ForceField.perceiveMMFF94AromaticRings(p, True)
ForceField.assignMMFF94AtomTypes(p, True, True)
ForceField.assignMMFF94BondTypeIndices(p, True, True)
ForceField.calcMMFF94AtomCharges(p, True, True)


# define positive and negative probe
positiveInteractionMatrix = np.zeros((len(randomSurfacePoints), p.numAtoms))  # interaction energy of positive probe at point x, y, z with all atoms
negativeInteractionMatrix = np.zeros((len(randomSurfacePoints), p.numAtoms))  # interaction energy of negative probe at point x, y, z with all atoms

scale_fact = 1.0
de_const = 1.0
dist_expo = 1.0

# calculate electrostatic potentials between probe at position x, y, z and atom n
for i in range(len(randomSurfacePoints)):
    point = randomSurfacePoints[i]
    for atomIndex in range(p.numAtoms):
        a = p.getAtom(atomIndex)

        partialCharge = ForceField.MMFF94Charge(a)
        negativeEnergy = ForceField.calcMMFF94ElectrostaticEnergy(coords[atomIndex], point, -1.0, partialCharge, scale_fact, de_const, dist_expo)
        positiveEnergy = ForceField.calcMMFF94ElectrostaticEnergy(coords[atomIndex], point, 1.0, partialCharge, scale_fact, de_const, dist_expo)

        negativeInteractionMatrix[i, atomIndex] = negativeEnergy
        positiveInteractionMatrix[i, atomIndex] = positiveEnergy
