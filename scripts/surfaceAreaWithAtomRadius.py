import numpy as np
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import CDPL.Pharm as Pharm
import CDPL.Base as Base
import CDPL.Shape as Shape
import sys
sys.path.append('../utils')
from Protein import Protein
from MoleculeTools import clean, sanitize_mol

def getShape(mol, scaleFactor=1):
    shape = Shape.GaussianShape()
    shapeFunc = Shape.GaussianShapeFunction()
    Shape.generateGaussianShape(mol, shape, inc_h=True)
    for e in shape:
        e.setRadius(e.getRadius()*scaleFactor)
    shapeFunc.setMaxOrder(6)
    shapeFunc.setShape(shape)
    return shape, shapeFunc

def getShapeWithIncreasedRadius(mol, increase=0.5):
    shape = Shape.GaussianShape()
    shapeFunc = Shape.GaussianShapeFunction()
    Shape.generateGaussianShape(mol, shape, inc_h=True)
    for e in shape:
        e.setRadius(e.getRadius()+increase)
    shapeFunc.setMaxOrder(6)
    shapeFunc.setShape(shape)
    return shape, shapeFunc


path = '../Data/benchmark5.5/structures/'
p = Protein()
p.fromFile('{}1A2K_l_b.pdb'.format(path))
# remove ligands and other crystalization artifacts
p.removeLigands()
sanitized = sanitize_mol(p, makeHydrogenComplete=True)
Pharm.prepareForPharmacophoreGeneration(p)
Chem.generateHydrogen3DCoordinates(p, True)


for i in range(10):
    scale = 1+i/10
    shape, shapeFunc = getShape(p, scaleFactor=scale)
    print(scale, shapeFunc.surfaceArea)