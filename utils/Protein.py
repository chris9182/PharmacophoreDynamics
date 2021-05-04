import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Shape as Shape
import CDPL.Math as Math
from typing import List
from scipy.spatial.transform import Rotation

from PharmacophoreTools import get_pharmacophore
from ProteinTools import THREE_LETTER_AMINO_ACID_CODES
from MathTools import *


class Protein(Chem.BasicMolecule):

    def __init__(self, structure: Chem.BasicMolecule = None):
        self.shape: Shape.GaussianShape = Shape.GaussianShape()
        self.shapeFunc: Shape.GaussianShapeFunction = Shape.GaussianShapeFunction()
        self.coordinates: Math.Vector3DArray = Math.Vector3DArray()
        self.ligands: List[Chem.BasicMolecule] = []

        super(Protein, self).__init__()
        if structure:
            self.assign(structure)

    def fromFile(self, path: str) -> Chem.BasicMolecule:
        from ProteinTools import readPDBFromFile

        self.assign(readPDBFromFile(path))
        return self

    def fromString(self, string: str) -> Chem.BasicMolecule:
        from ProteinTools import readPDBFromString

        self.assign(readPDBFromString(string))
        return self

    def fromRSCB(self, pdbCode: str) -> Chem.BasicMolecule:
        from ProteinTools import readPDBFromRSCB

        self.assign(readPDBFromRSCB(pdbCode))
        return self

    def separateLigandFromProtein(self, keep: bool = True, removeWater: bool = True) -> None:
        """
        Removes all entities from the protein which are not an amino acid --> usually a ligand.
        Be careful with peptides and covalently bound ligands! These will face a different behaviour and might cause
        unexpected results!
        :param keep: Whether to keep the removed ligands stored in the ligand property of just remove them.
        :param removeWater: If true, removes the water molecules
        :return:
        """
        from ProteinTools import getMoleculeFromAtom

        atomsToRemoveFromProtein = []
        for atom in self.atoms:
            index = self.getAtomIndex(atom)
            if index in atomsToRemoveFromProtein:
                continue

            ligandCode = Biomol.getResidueCode(atom)

            if ligandCode == "HOH":
                if removeWater:
                    atomsToRemoveFromProtein.append(index)
                    continue
                else:
                    raise NotImplementedError

            if ligandCode not in THREE_LETTER_AMINO_ACID_CODES:
                ligand, atomsToRemove = getMoleculeFromAtom(atom, self)
                print(atomsToRemove)
                atomsToRemoveFromProtein.extend(atomsToRemove)
                if keep:
                    self.ligands.append(ligand)

        atomsToRemoveFromProtein = list(set(atomsToRemoveFromProtein))
        atomsToRemoveFromProtein.sort()
        for i in reversed(atomsToRemoveFromProtein):
            self.removeAtom(i)

    def getGaussianShape(self, addFeatures: bool = False, maxOrder: int = 1) -> Shape.GaussianShape:
        Shape.generateGaussianShape(self, self.shape)
        if addFeatures:
            ph4 = get_pharmacophore(self)
            Shape.generateGaussianShape(ph4, self.shape, True)  # add pharmacophore shape to molecule shape
        self.shapeFunc.setMaxOrder(maxOrder)
        self.shapeFunc.setShape(self.shape)

        return self.shape

    def makeRandomRotation(self, inplace: bool = True) -> Math.Vector3DArray:
        # TODO: maybe add boundaries for randomness
        rotMatrix = Math.Matrix3D()
        rotMatrix.assign(Rotation.random().as_matrix())
        rotatedCoords = rotate3DObject(self.getCoordinates(), rotMatrix)

        if inplace:
            Chem.set3DCoordinates(rotatedCoords)

        return rotatedCoords

    def makeRandomTranslation(self, inplace: bool = True, scalingFactor: float = 10) -> Math.Vector3DArray:
        """
        
        :param inplace: 
        :param scalingFactor: Scales the randomly retrieved direction by this factor
        :return: 
        """
        direction = Math.Vector3D()
        direction.assign(np.random.rand(3)*scalingFactor)
        translatedCoords = translate3DObject(self.getCoordinates(), direction)

        if inplace:
            Chem.set3DCoordinates(translatedCoords)

        return translatedCoords

    def getSurfaceExposedAtoms(self) -> List[Chem.BasicAtom]:  # TODO
        """
        Get a list of CDPL Molecules located on the protein surface.
        :return:
        """
        raise NotImplementedError

    def generateSurfacePoints(self, scaleFactor: float = 1, density: float = 1) -> Math.Vector3DArray:
        """
        Get an array of coordinates corresponding to surface exposed atoms in the protein. Points are placed on the vdw
        surface of the atoms.
        :param scaleFactor: Scales the VDW radius. Points are set at a distance proportional to this factor.
        :param density: Determines the density of points representing the atom surface per square angstrom.
        :return:
        """
        from pyvdwsurface import vdwsurface

        atomTypes = [Chem.getSymbol(a) for a in self.atoms]
        points = vdwsurface(self.getCoordinates(), atomTypes, scale_factor=scaleFactor, double_density=density)
        surfacePointCoordinates = Math.Vector3DArray()
        surfacePointCoordinates.assign(points)
        return surfacePointCoordinates

    def getCoordinates(self) -> Math.Vector3DArray:
        if self.coordinates.isEmpty():
            Chem.get3DCoordinates(self, self.coordinates)
        return self.coordinates
