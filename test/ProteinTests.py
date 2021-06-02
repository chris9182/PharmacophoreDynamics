import unittest
import os
from numpy.random import random_sample
import CDPL.Chem as Chem

from src.utils.Protein import Protein
from test.TestConstants import *


class ProteinTests(unittest.TestCase):


    def setUp(self) -> None:
        self.protein = loadCleanedProtein()

    def loadPDBTest(self):
        p = loadProtein()
        self.assertEqual(type(p), Chem.BasicMolecule)
        print(p.numAtoms)
        self.assertEqual(p.numAtoms, NUMBER_OF_ATOMS)

    def cleanLigandTest(self):
        p = loadProtein()

        # assert ligand and other artefacts are included
        numAtomsBefore = p.numAtoms
        self.assertEqual(numAtomsBefore, NUMBER_OF_ATOMS)
        numberOfEntitiesBefore = Chem.getComponentCount(p)
        self.assertGreater(numberOfEntitiesBefore, 1)

        p.removeLigands()  # clean protein

        # assert successful cleaning
        self.assertLess(p.numAtoms, numAtomsBefore)
        self.assertEqual(p.numAtoms, NUMBER_OF_ATOMS_WITHOUT_LIGAND)
        self.assertEqual(Chem.getComponentCount(p), 1)
        self.assertLess(Chem.getComponentCount(p), numberOfEntitiesBefore)

    def addHydrogensTest(self):
        p = loadProtein()
        p.removeLigands()

        numAtomsBefore = p.numAtoms
        self.assertEqual(numAtomsBefore, NUMBER_OF_ATOMS_WITHOUT_LIGAND)

        p.addHydrogens()
        self.assertGreater(p.numAtoms, numAtomsBefore)
        self.assertEqual(p.numAtoms, NUMBER_OF_ATOMS_WITH_HYDROGENS_CLEANED)

    @unittest.skip('not implemented')
    def writePDBTest(self):
        raise NotImplementedError


def loadProtein() -> Protein:
    p = Protein()
    return p.fromFile('{}{}_l_b.pdb'.format(DATA_FOLDER, PDB_CODE))


def loadCleanedProtein() -> Protein:
    p = loadProtein()
    p.prepare()
    return p