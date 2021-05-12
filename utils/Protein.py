import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Shape as Shape
import CDPL.Math as Math
from typing import List
from scipy.spatial.transform import Rotation

from PharmacophoreTools import getPharmacophore
from ProteinTools import THREE_LETTER_AMINO_ACID_CODES
from MathTools import *


class Protein(Chem.BasicMolecule):

    def __init__(self, structure: Chem.BasicMolecule = None):
        self.shape: Shape.GaussianShape = Shape.GaussianShape()
        self.shapeFunc: Shape.GaussianShapeFunction = Shape.GaussianShapeFunction()
        self.coordinates: Math.Vector3DArray = Math.Vector3DArray()
        self.ligands: List[Chem.BasicMolecule] = []
        self.surfaceAtoms: Chem.BasicMolecule = Chem.BasicMolecule()

        super(Protein, self).__init__()
        if structure:
            from MoleculeTools import sanitize_mol
            self.assign(structure)
            # sanitize_mol(self, makeHydrogenComplete=True)

    def prepare(self, removeLigands=True):
        from MoleculeTools import sanitize_mol

        sanitize_mol(self, makeHydrogenComplete=True)
        Chem.generateHydrogen3DCoordinates(self, True)

        if removeLigands:
            self.removeLigands()

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

    def toFile(self, path: str) -> str:
        from ProteinTools import writePDB

        writePDB(path, self)
        return path

    def removeLigands(self, keep: bool = False, removeWater: bool = True) -> None:
        """
        Removes all entities from the protein which are not an amino acid --> usually a ligand.
        Be careful with peptides and covalently bound ligands! These will face a different behaviour and might cause
        unexpected results!
        :param keep: Whether to keep the removed ligands stored in the ligand property or simply remove them.
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

    def getGaussianShape(self, addFeatures: bool = False, maxOrder: int = 4) -> Shape.GaussianShape:
        self.shape.clear()
        Shape.generateGaussianShape(self, self.shape)
        if addFeatures:
            ph4 = getPharmacophore(self)
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
            Chem.set3DCoordinates(self, rotatedCoords)

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
            Chem.set3DCoordinates(self, translatedCoords)

        return translatedCoords

    def getSurfaceExposedAtoms(self, copy=True) -> Chem.BasicMolecule:
        """
        Get a list of CDPL Molecules located on the protein surface.
        :param copy: Whether to return a copy of the surface atoms or the surface atom object itself.
        :return:
        """
        from ProteinTools import getSurfaceAtoms

        self.surfaceAtoms.assign(getSurfaceAtoms(self))

        if copy:
            surfaceAtoms = Chem.BasicMolecule()
            surfaceAtoms.assign(self.surfaceAtoms)
            return  surfaceAtoms

        return self.surfaceAtoms

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
