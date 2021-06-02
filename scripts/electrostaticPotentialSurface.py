import CDPL.Math as Math
import CDPL.Chem as Chem
import CDPL.ForceField as ForceField
import CDPL.Biomol as Biomol
import CDPL.Util as Util
import numpy as np

from src.utils.Protein import Protein

# load protein
p = Protein()
p.fromFile('../data/../data/benchmark5.5/structures/1A2K_l_b.pdb')
p.prepare()  # calculates H-Atoms!
coords = p.getCoordinates().toArray(False)

randomSurfacePoints = np.random.randint(0, 100, size=(100, 3))

# prepare protein for force field calculations
sssr = ForceField.MMFF94AromaticSSSRSubset(p)
ForceField.setMMFF94AromaticRings(p, sssr)
ForceField.assignMMFF94AtomTypes(p, True, True)
ForceField.assignMMFF94BondTypeIndices(p, True, True)
ForceField.calcMMFF94AtomCharges(p, True, True)

# calculate partial charges
chargeCalculator = ForceField.MMFF94ChargeCalculator()
partialCharges = Util.DArray()
chargeCalculator.calculate(p, partialCharges, True)
partialCharges = np.array(partialCharges)

# define positive and negative probe
probeCharges = [-1, 1]
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
        negativeElectrosstaticInteraction = ForceField.MMFF94ElectrostaticInteraction(0,  # probe index
                                                                                      atomIndex,
                                                                                      probeCharges[0],  # charge of probe
                                                                                      partialCharges[atomIndex],
                                                                                      scale_fact,  # no idea what this does
                                                                                      de_const,   # no idea what this does
                                                                                      dist_expo  # no idea what this does
                                                                                      )
        positiveElectrosstaticInteraction = ForceField.MMFF94ElectrostaticInteraction(1,  # probe index
                                                                                      atomIndex,
                                                                                      probeCharges[1],  # charge of probe
                                                                                      partialCharges[atomIndex],
                                                                                      scale_fact,  # no idea what this does
                                                                                      de_const,  # no idea what this does
                                                                                      dist_expo  # no idea what this does
                                                                                      )

        atomCoordinates = coords[atomIndex]
        coordinatePair = np.array([coords[atomIndex], point])

        negativeEnergy = ForceField.calcMMFF94ElectrostaticEnergy(negativeElectrosstaticInteraction, coordinatePair)
        positiveEnergy = ForceField.calcMMFF94ElectrostaticEnergy(positiveElectrosstaticInteraction, coordinatePair)
        negativeInteractionMatrix[i, atomIndex] = negativeEnergy
        positiveInteractionMatrix[i, atomIndex] = positiveEnergy
