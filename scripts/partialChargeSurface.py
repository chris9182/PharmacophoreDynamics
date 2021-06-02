import CDPL.Math as Math
import CDPL.Chem as Chem
import CDPL.ForceField as ForceField
import CDPL.Biomol as Biomol
import CDPL.Util as Util
import numpy as np
from sklearn.preprocessing import MinMaxScaler

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

# get partial charge for surface point as weighted average of atoms
surfaceCharges = []
for pointCoords in randomSurfacePoints:
    distances = np.linalg.norm(coords-pointCoords, ord=2, axis=1)
    weights = 1/distances
    surfacePointPartialCharge = np.average(partialCharges, weights=weights)
    surfaceCharges.append(surfacePointPartialCharge)

# scale partial charges to [-1, 1]
scaler = MinMaxScaler(feature_range=(-1, 1))
scaledCharges = scaler.fit_transform(np.array(surfaceCharges).reshape(-1, 1))
print(scaledCharges)
