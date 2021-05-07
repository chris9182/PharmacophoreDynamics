import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Base as Base


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


SINGLE_LETTER_AMINO_ACID_CODES = (  # same order as three letter codes
    'A',  # alanine
    'R',  # arginine
    'N',  # asparagine
    'D',  # aspartic acid
    'B',  # aspartic acid or asparagine
    'C',  # cysteine
    'E',  # glutamic acid
    'Q',  # glutamine
    'Z',  # glutamic acid or glutamine
    'G',  # glycine
    'H',  # histidin
    'I',  # isoleucine
    'L',  # leucine
    'K',  # lysine
    'M',  # methionine
    'F',  # phenylalanine
    'P',  # proline
    'S',  # serine
    'T',  # threonine
    'W',  # tryptophan
    'Y',  # tyrosine
    'V',  # valine
)


def readPDBFromStream(stream: Base.IOStream):
    from Protein import Protein
    from MoleculeTools import sanitize_mol

    r = Biomol.PDBMoleculeReader(stream)
    mol = Chem.BasicMolecule()
    r.read(mol)
    sanitize_mol(mol, makeHydrogenComplete=True)
    return Protein(mol)


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


def writePDB(path: str) -> None:
    w = Biomol.FilePDBMolecularGraphWriter(path)
    w.write()
    w.close()


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
