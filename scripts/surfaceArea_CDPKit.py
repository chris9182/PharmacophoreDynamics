import CDPL.Chem as Chem
import CDPL.Shape as Shape
import CDPL.Math as Math
import CDPL.Biomol as Biomol
import CDPL.Base as Base

from typing import List


THREE_LETTER_AMINO_ACID_CODES = (
    'ALA',
    'ARG',
    'ASN',
    'ASP',
    'ASX',  # asparagine or aspartic acid
    'CYS',
    'GLU',
    'GLN',
    'GLX',  # glutamine or glutamic acid
    'GLY',
    'HIS',
    'ILE',
    'LEU',
    'LYS',
    'MET',
    'PHE',
    'PRO',
    'SER',
    'THR',
    'TRP',
    'TYR',
    'VAL'
)


def getGaussianShapeOfMolecule(
        mol: Chem.BasicMolecule,
        addFeatures: bool = False,
        maxOrder: int = 1) -> (Shape.GaussianShape, Shape.GaussianShapeFunction):
    shape = Shape.GaussianShape()
    shapeFunc = Shape.GaussianShapeFunction()
    Shape.generateGaussianShape(mol, shape)
    # if addFeatures:
    #     ph4 = getPharmacophore(mol)
    #     Shape.generateGaussianShape(ph4, shape, True)  # add pharmacophore shape to molecule shape
    shapeFunc.setMaxOrder(maxOrder)
    shapeFunc.setShape(shape)

    return shape, shapeFunc


def readPDBFromStream(stream: Base.IOStream):
    r = Biomol.PDBMoleculeReader(stream)
    mol = Chem.BasicMolecule()
    r.read(mol)
    return Protein(mol)


def readPDBFromFile(path: str) -> Chem.BasicMolecule:
    s = Base.FileIOStream(path)
    protein = readPDBFromStream(s)
    return protein


def getMoleculeFromAtom(atom: Chem.BasicAtom, protein: Chem.BasicMolecule) -> (Chem.BasicMolecule, list):
    """
    Given an atom and a protein structure, find the ligand the atom corresponds to.
    Traverses the molecule by its bonds until no longer any atoms are attached. All atoms and bonds are assigned to a
    new molecule, which is being returned.
    :param atom:
    :param protein:
    :return: The found ligand as well as the atom indices of the ligand in the parent molecule.
    """
    ligand = Chem.Fragment()
    neighbors = set()  # atoms not being added already
    neighborsAdded = set()  # keep track of added atoms to not process twice
    atomsToRemove = []

    neighbors.add(atom)
    while len(neighbors) > 0:
        n = neighbors.pop()
        neighborsAdded.add(n)
        ligand.addAtom(n)
        atomsToRemove.append(protein.getAtomIndex(n))

        # get all the neighbor atoms
        for i, b in enumerate(n.bonds):
            for a in b.atoms:
                if a != n:
                    if a not in neighbors and a not in neighborsAdded:  # new atom
                        neighbors.add(a)

                    ligand.addBond(b)  # ignored if already exists

    Chem.perceiveComponents(ligand, True)
    mol = Chem.BasicMolecule()
    mol.assign(ligand)
    return mol, atomsToRemove



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
        self.assign(readPDBFromFile(path))
        return self

    # def fromString(self, string: str) -> Chem.BasicMolecule:
    #     from ProteinTools import readPDBFromString
    #
    #     self.assign(readPDBFromString(string))
    #     return self
    #
    # def fromRSCB(self, pdbCode: str) -> Chem.BasicMolecule:
    #     from ProteinTools import readPDBFromRSCB
    #
    #     self.assign(readPDBFromRSCB(pdbCode))
    #     return self

    def separateLigandFromProtein(self, keep: bool = True, removeWater: bool = True) -> None:
        """
        Removes all entities from the protein which are not an amino acid --> usually a ligand.
        Be careful with peptides and covalently bound ligands! These will face a different behaviour and might cause
        unexpected results!
        :param keep: Whether to keep the removed ligands stored in the ligand property of just remove them.
        :param removeWater: If true, removes the water molecules
        :return:
        """
        # from ProteinTools import getMoleculeFromAtom

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
                # print(atomsToRemove)
                atomsToRemoveFromProtein.extend(atomsToRemove)
                if keep:
                    self.ligands.append(ligand)

        atomsToRemoveFromProtein = list(set(atomsToRemoveFromProtein))
        atomsToRemoveFromProtein.sort()
        for i in reversed(atomsToRemoveFromProtein):
            self.removeAtom(i)

    def getGaussianShape(self, addFeatures: bool = False, maxOrder: int = 1) -> Shape.GaussianShape:
        Shape.generateGaussianShape(self, self.shape)
        # if addFeatures:
        #     ph4 = getPharmacophore(self)
        #     Shape.generateGaussianShape(ph4, self.shape, True)  # add pharmacophore shape to molecule shape
        self.shapeFunc.setMaxOrder(maxOrder)
        self.shapeFunc.setShape(self.shape)

        return self.shape

    # def makeRandomRotation(self, inplace: bool = True) -> Math.Vector3DArray:
    #     # TODO: maybe add boundaries for randomness
    #     rotMatrix = Math.Matrix3D()
    #     rotMatrix.assign(Rotation.random().as_matrix())
    #     rotatedCoords = rotate3DObject(self.getCoordinates(), rotMatrix)
    #
    #     if inplace:
    #         Chem.set3DCoordinates(rotatedCoords)
    #
    #     return rotatedCoords
    #
    # def makeRandomTranslation(self, inplace: bool = True, scalingFactor: float = 10) -> Math.Vector3DArray:
    #     """
    #
    #     :param inplace:
    #     :param scalingFactor: Scales the randomly retrieved direction by this factor
    #     :return:
    #     """
    #     direction = Math.Vector3D()
    #     direction.assign(np.random.rand(3) * scalingFactor)
    #     translatedCoords = translate3DObject(self.getCoordinates(), direction)
    #
    #     if inplace:
    #         Chem.set3DCoordinates(translatedCoords)
    #
    #     return translatedCoords

    def getSurfaceExposedAtoms(self) -> List[Chem.BasicAtom]:  # TODO
        """
        Get a list of CDPL Molecules located on the protein surface.
        :return:
        """
        raise NotImplementedError

    # def generateSurfacePoints(self, scaleFactor: float = 1, density: float = 1) -> Math.Vector3DArray:
    #     """
    #     Get an array of coordinates corresponding to surface exposed atoms in the protein. Points are placed on the vdw
    #     surface of the atoms.
    #     :param scaleFactor: Scales the VDW radius. Points are set at a distance proportional to this factor.
    #     :param density: Determines the density of points representing the atom surface per square angstrom.
    #     :return:
    #     """
    #     from pyvdwsurface import vdwsurface
    #
    #     atomTypes = [Chem.getSymbol(a) for a in self.atoms]
    #     points = vdwsurface(self.getCoordinates(), atomTypes, scale_factor=scaleFactor, double_density=density)
    #     surfacePointCoordinates = Math.Vector3DArray()
    #     surfacePointCoordinates.assign(points)
    #     return surfacePointCoordinates

    def getCoordinates(self) -> Math.Vector3DArray:
        if self.coordinates.isEmpty():
            Chem.get3DCoordinates(self, self.coordinates)
        return self.coordinates


if __name__ == '__main__':
    # load protein
    protein = Protein()
    protein.fromFile('1A2K_l_b.pdb')

    # remove ligand and water, just to be sure
    protein.separateLigandFromProtein()

    # calculate surface area for a single carbon atom
    protein.getGaussianShape()
    shapeFunc = protein.shapeFunc
    carbonProteinIndex = None
    for a in protein.atoms:
        if Chem.getType(a) == 6:
            carbonProteinIndex = protein.getAtomIndex(a)
            break
    surfAreaCarbonProtein = shapeFunc.calcSurfaceArea(carbonProteinIndex)  # calculate the surface area contribution here?

    # create a simple carbon molecule with coordinates
    carbon = Chem.BasicMolecule()
    cAtom = carbon.addAtom()
    Chem.setType(cAtom, 6)
    coords = Math.Vector3D()
    coords.assign([1, 2, 3])
    Chem.set3DCoordinates(cAtom, coords)

    # calculate shape and surface area of carbon molecule
    carbonShape, carbonShapeFunc = getGaussianShapeOfMolecule(carbon)
    surfAreaCarbon = carbonShapeFunc.surfaceArea

    # assert that contribution of carbon atom in protein is in fact the entire surface area of a single carbon
    assert surfAreaCarbon == surfAreaCarbonProtein

    # What am I missing here?
    # Summing the surface area of all atoms in the protein yields the surface area of the protein. I find it hard to
    # believe that all the atoms are surface atoms. Am I missing the VDW radius somehow and all atoms are in fact
    # just points right now?
    # for example:
    # shapeFunc.distCutoff -> 0
    # Already tried: shapeFunc.setDistanceCutoff(1) --> still same surface area and the same contribution.


    # calculate the surface area of the protein
    surfaceAreaProtein = protein.shapeFunc.surfaceArea

    # now do it again, but sum the contributions myself
    totalSA = 0
    indivSurfaceAreas = []
    for a in protein.atoms:
        sa = protein.shapeFunc.calcSurfaceArea(protein.getAtomIndex(a))
        indivSurfaceAreas.append(sa)
        assert sa > 0  # assume all atoms are on the surface --> contribution is larger than 0. In fact, they are all ~ 30 Angstrom²
        totalSA += sa

    assert surfaceAreaProtein == totalSA  # surface area is sum of ALL atoms, even if not on surface

    # print(indivSurfaceAreas)
