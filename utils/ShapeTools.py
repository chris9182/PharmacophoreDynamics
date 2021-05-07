import CDPL.Shape as Shape
import CDPL.Math as Math
import CDPL.Chem as Chem
from PharmacophoreTools import getPharmacophore


def getGaussianShapeOfMolecule(
        mol: Chem.BasicMolecule,
        addFeatures: bool = False,
        maxOrder: int = 4) -> (Shape.GaussianShape, Shape.GaussianShapeFunction):
    shape = Shape.GaussianShape()
    shapeFunc = Shape.GaussianShapeFunction()
    Shape.generateGaussianShape(mol, shape)
    if addFeatures:
        ph4 = getPharmacophore(mol)
        Shape.generateGaussianShape(ph4, shape, True)  # add pharmacophore shape to molecule shape
    shapeFunc.setMaxOrder(maxOrder)
    shapeFunc.setShape(shape)

    return shape, shapeFunc
