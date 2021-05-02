import CDPL.Biomol as Biomol
import CDPL.Chem as Chem
import CDPL.Shape as Shape
from PharmacophoreTools import get_pharmacophore
from ProteinTools import THREE_LETTER_AMINO_ACID_CODES


class Protein(Chem.BasicMolecule):

    def __init__(self, structure: Chem.BasicMolecule = None):
        self.shape = None
        self.shapeFunc = None
        self.ligands = []

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

    def separateLigandFromProtein(self, keep=True, removeWater=True):
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

    def getGaussianShape(self, addFeatures: bool =False, maxOrder: int =1, recalculate: bool =True) -> Shape.GaussianShape:
        if not self.shape or recalculate:
            self.shape = Shape.GaussianShape()
            self.shapeFunc = Shape.GaussianShapeFunction()
            Shape.generateGaussianShape(self, self.shape)
            if addFeatures:
                ph4 = get_pharmacophore(self)
                Shape.generateGaussianShape(ph4, self.shape, True)  # add pharmacophore shape to molecule shape
            self.shapeFunc.setMaxOrder(maxOrder)
            self.shapeFunc.setShape(self.shape)

        return self.shape
