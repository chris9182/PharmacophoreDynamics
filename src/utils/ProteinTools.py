import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Base as Base


def readPDBFromStream(stream: Base.IOStream):
    from src.utils.Protein import Protein
    # from src.utils.MoleculeTools import sanitize_mol

    r = Biomol.PDBMoleculeReader(stream)
    p = Protein()
    r.read(p)
    # sanitize_mol(mol, makeHydrogenComplete=True)
    return p


def readPDBFromFile(path: str) -> Chem.BasicMolecule:
    s = Base.FileIOStream(path)
    protein = readPDBFromStream(s)
    return protein


def readPDBFromString(string: str) -> Chem.BasicMolecule:
    s = Base.StringIOStream(string)
    protein = readPDBFromStream(s)
    return protein


def readPDBFromRSCB(pdbCode: str) -> Chem.BasicMolecule:
    from requests import request

    r = request('GET', 'https://files.rcsb.org/download/{pdb}.pdb'.format(pdb=pdbCode))
    if r.ok:
        return readPDBFromString(r.text)
    else:
        return None


def writePDB(path: str, protein: Chem.BasicMolecule) -> None:
    Chem.makeHydrogenDeplete(protein)
    w = Biomol.FilePDBMolecularGraphWriter(path)
    w.write(protein)
    w.close()


# def prepareProtein(protein, removeLigands=True, removeWater=True):
#     from src.utils.MoleculeTools import sanitize_mol
#
#     sanitize_mol(self, makeHydrogenComplete=True)
#     Chem.generateHydrogen3DCoordinates(self, True)
#
#     if removeLigands:
#         self.removeLigands(removeWater=removeWater)
#
#     return protein


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


def getSurfaceAtoms(mol):
    surfaceATomExtractor = Chem.SurfaceAtomExtractor()
    f = Chem.Fragment()
    surfaceATomExtractor.extract(mol, mol, f)
    surfaceAtoms = Chem.BasicMolecule()
    surfaceAtoms.assign(f)
    return surfaceAtoms
