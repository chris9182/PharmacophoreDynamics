import sys
sys.path.append("/data/shared/software/CDPKit/Python")
import CDPL.Chem as Chem
import CDPL.Base as Base
import CDPL.Util as Util
import CDPL.Math as Math
import CDPL.Pharm as Pharm
import numpy as np
from PharmacophoreTools import getPharmacophore
from collections import Iterable


ALLOWED_ATOMS = [1, 6, 7, 8, 9, 15, 16, 17, 35, 53]
SALT_METALS = [3, 11, 12, 19, 20]


class SDFReader:

    """
    Handles reading molecules from SDF files and takes care of necessary preparation steps. 
    """

    def __init__(self, path: str, multiconf: bool=True, nr_mols: int=-1, properties: list=None):
        """

        :param path:
        :param multiconf:
        :param nr_mols:
        :param properties: List of properties to be read from SDF file. Returns None if not found. Internally adds ' <'
        and '>' before and after the property names in order to comply with the toolkit.
        """
        self.r = Chem.FileSDFMoleculeReader(path)
        self.multiconf = multiconf
        self.nr_mols = nr_mols
        self.nr_samples = None
        self.properties = set([" <"+p+">" for p in properties]) if properties is not None else None

        Chem.setMultiConfImportParameter(self.r, multiconf)
        if multiconf:
            Chem.setMultiConfInputProcessorParameter(self.r, Chem.DefaultMultiConfMoleculeInputProcessor(False,
                                                                                                        Chem.AtomPropertyFlag.TYPE | Chem.AtomPropertyFlag.ISOTOPE | Chem.AtomPropertyFlag.FORMAL_CHARGE,
                                                                                                        Chem.BondPropertyFlag.ORDER))
        self._gen = iter(self)

    def __len__(self):
        if self.nr_samples is None:
            self.nr_samples = self.r.getNumRecords()
            return self.nr_samples
        return self.nr_samples

    def __iter__(self):
        if self.properties is None:
            i = 0
            while True:
                mol = Chem.BasicMolecule()
                try:
                    if self.r.read(mol):
                        yield sanitize_mol(mol)
                    else:
                        break
                except IOError:
                    yield None

                i += 1
                if i == self.nr_mols:
                    break
        else:
            i = 0
            while True:
                mol = Chem.BasicMolecule()
                try:
                    if self.r.read(mol):
                        read_properties = self._extract_properties_from_mol(mol)
                        yield sanitize_mol(mol), read_properties
                    else:
                        break
                except IOError:
                    yield None

                i += 1
                if i == self.nr_mols:
                    break

    def __call__(self):
        return iter(self)

    def _extract_properties_from_mol(self, mol):
        read_properties = {}
        data = Chem.getStructureData(mol)
        for element in data:
            if element.header in self.properties:
                read_properties[element.header[2:-1]] = element.data
        return read_properties

    def read_all(self):
        """
        Reads all the molecules from the SDF file with set properties
        :return:
        """
        mols = {}
        for i, mol in enumerate(self):
            name = Chem.getName(mol)
            if len(name) == 0:  # no name set
                name = str(i)
            mols[name] = mol
        return mols

    def read(self):
        try:
            return next(self._gen)
        except IOError:
            return None


def remove_metal_salts(mol: Chem.BasicMolecule) -> Chem.BasicMolecule:
    to_remove = []
    for atom in mol.atoms:
        if Chem.getType(atom) not in SALT_METALS:
            continue
        else:
            to_remove.append(mol.getAtomIndex(atom))
    to_remove.sort()
    to_remove.reverse()
    for index in to_remove:
        mol.removeAtom(index)
    return mol


def is_metal(mol: Chem.BasicMolecule) -> bool:
    """
    Indicate if the compound contains a metal
    """
    for atom in mol.atoms:
        if Chem.getType(atom) in ALLOWED_ATOMS:
            continue
        else:
            return True
    return False


def remove_components(mol: Chem.BasicMolecule) -> Chem.BasicMolecule:
    components = Chem.getComponents(mol)
    largest_component = None  # set default number of atoms and index
    for comp in components:
        if largest_component is None:
            largest_component = comp
        elif comp.numAtoms > largest_component.numAtoms:
            largest_component = comp
    new_mol = Chem.BasicMolecule()
    new_mol.assign(largest_component)
    if Chem.hasStructureData(mol):
        Chem.setStructureData(new_mol, Chem.getStructureData(mol))
    return new_mol


def is_inorganic(mol: Chem.BasicMolecule) -> bool:
    for atom in mol.atoms:
        if Chem.getType(atom) != 6:
            continue
        else:
            return False
    return True


def neutralise(mol: Chem.BasicMolecule) -> Chem.BasicMolecule:
    to_remove = []
    for atom in mol.atoms:
        if Chem.getFormalCharge(atom) != 0:
            form_charge = Chem.getFormalCharge(atom)

            if form_charge != 0:
                for nbr_atom in atom.atoms:
                    if Chem.getFormalCharge(nbr_atom) != 0:
                        form_charge = 0
                        break  # it's fine if neighbor is charged too -> we assume it's the opposite charge

            if form_charge != 0:
                if form_charge > 0:
                    form_charge -= Chem.getImplicitHydrogenCount(atom)

                    if form_charge < 0:  # if charge is negative we set to zero and calculate the number of hydrogens later on
                        form_charge = 0

                    for nbr_atom in atom.atoms:
                        if form_charge == 0:
                            break

                        if Chem.getType(nbr_atom) == Chem.AtomType.H:
                            to_remove.append(mol.getAtomIndex(nbr_atom))
                            form_charge -= 1

                    Chem.setFormalCharge(atom, form_charge)

                else:
                    Chem.setFormalCharge(atom, 0)

    if len(to_remove) > 0:
        to_remove.sort()
        to_remove.reverse()
        for index in to_remove:
            mol.removeAtom(index)

        for atom in mol.atoms:
            Chem.setImplicitHydrogenCount(atom, Chem.calcImplicitHydrogenCount(atom, mol))
            Chem.setHybridizationState(atom, Chem.perceiveHybridizationState(atom, mol))
    return mol


def sanitize_mol(mol: Chem.BasicMolecule, makeHydrogenComplete=False) -> Chem.BasicMolecule:
    Chem.calcImplicitHydrogenCounts(mol, True)
    Chem.perceiveHybridizationStates(mol, True)
    Chem.perceiveComponents(mol, True)
    Chem.perceiveSSSR(mol, True)
    Chem.setRingFlags(mol, True)
    Chem.setAromaticityFlags(mol, True)
    if makeHydrogenComplete:
        Chem.makeHydrogenComplete(mol)
        Chem.calcImplicitHydrogenCounts(mol, True)
        Chem.generateHydrogen3DCoordinates(mol, True)
    return mol


def clean(mol: Chem.BasicMolecule) -> "bool or Chem.BasicMolecule":
    """
    Checks if the molecule is a metal or inorganic -> return False
    Removes salts [Na+, Mg2+, Ca2+, K+, Li+] and multiple components -> keep the largest component.
    Neutralise molecule by adding or removing protons as well as possible.
    :param mol:
    :return:
    """
    mol = sanitize_mol(mol)

    # clean molecule
    mol = remove_metal_salts(mol)
    if is_metal(mol):
        return False
    if Chem.getComponentCount(mol) > 1:
        mol = remove_components(mol)
    if is_inorganic(mol):
        return False
    mol = neutralise(mol)
    return sanitize_mol(mol)


def mol_from_smiles(smiles: str) -> Chem.BasicMolecule:
    ifs = Base.StringIOStream(smiles)
    mol = Chem.BasicMolecule()
    r = Chem.SMILESMoleculeReader(ifs)
    r.read(mol)
    return sanitize_mol(mol)


def mol_from_smi(path: str) -> Chem.BasicMolecule:
    smiles = np.loadtxt(path, delimiter=",", dtype=str, comments=None)
    molecules = []
    for smile in smiles:
        mol = mol_from_smiles(smile)
        molecules.append(mol)
    if len(molecules) > 1:
        return molecules
    return molecules[0]


def calculateECFP(mol, nIter=4, nBits=1021):
    """
    Calculate the ECFP fingerprint for a given molecule.
    :param mol:
    :param nIter:
    :param nBits:
    :return:
    """
    Chem.makeHydrogenComplete(mol)
    ecfpGen = Chem.CircularFingerprintGenerator()
    ecfpGen.setNumIterations(nIter)
    ecfpGen.setNumBits(nBits)
    bitv = Util.BitSet()
    ecfpGen.generate(mol, bitv)
    return bitv


def calculate_tanimoto(fp1: np.ndarray, fp2:np.ndarray) -> float:
    nom = np.sum(fp1 * fp2)
    denom = np.sum(fp1) + np.sum(fp2) - nom
    return nom / denom


# def smiles_to_vector(smiles, mapping, zero_padding=False):
#     """
#     Transform the SMILES string into a zero-padded vector containing numbers corresponding ot the SMILES characters
#     according to the additional given mapping.
#     :param smiles:
#     :param mapping:
#     :return:
#     """
#     # translate SMILES
#     lengths = []
#     translated_smiles = []
#     for i, smile in enumerate(smiles):
#         smiles_vector = []
#         previous_char = None
#         for char in smile:
#             if char == "l":
#                 char = previous_char + char
#                 smiles_vector[-1] = mapping[char]
#             elif char == "r":
#                 continue
#             elif char == "B":
#                 char = "Br"
#                 smiles_vector.append(mapping[char])
#             else:
#                 smiles_vector.append(mapping[char])
#             previous_char = char
#         lengths.append(len(smiles_vector))
#         translated_smiles.append(smiles_vector)
#
#     if not zero_padding:
#         return translated_smiles, lengths
#     # zero pad SMILES
#     zero_padded_vector = np.zeros((len(translated_smiles), max(lengths)), dtype=np.int16)
#     for i, smiles_vector in enumerate(translated_smiles):
#         l = len(smiles_vector)
#         zero_padded_vector[i, :l] = smiles_vector
#     return zero_padded_vector, np.array(lengths, dtype=np.int16)


def get_span(mol: Chem.BasicMolecule) -> int:
    Chem.calcTopologicalDistanceMatrix(mol, False)
    return Chem.calcTopologicalDiameter(mol)


def check_chemical_validity(mol):
    raise NotImplementedError


def get_centroid(coordinates: np.ndarray) -> (float, float, float):
    length = coordinates.shape[0]
    sum_x = np.sum(coordinates[:, 0])
    sum_y = np.sum(coordinates[:, 1])
    sum_z = np.sum(coordinates[:, 2])
    return sum_x / length, sum_y / length, sum_z / length


def center_mol(mol):
    coords = Math.Vector3DArray()
    Chem.get3DCoordinates(mol, coords)
    np_coords = np.array(coords)
    centroid = get_centroid(np_coords)
    centered = np_coords - centroid

    # set coordinates coordinate object
    for i, row in enumerate(coords):
        for j, column in enumerate(row):
            row[j] = centered[i, j]

    # set coordinates to molecule
    Chem.set3DCoordinates(mol, coords)
    return mol


def translate_mol_to_coords(mol, coords):
    mol_coords = Math.Vector3DArray()
    Chem.get3DCoordinates(mol, mol_coords)
    for i, row in enumerate(mol_coords):
        for j, column in enumerate(row):
            row[j] = column - coords[j]  # shift to desired coordinates

    Chem.set3DCoordinates(mol, mol_coords)
    return mol


def mol_to_smiles(mol, kekulized=False, canonical=True, atom_stereo=True, hydrogen_deplete=True, bond_stereo=False):
    stream = Base.StringIOStream()
    w = Chem.SMILESMolecularGraphWriter(stream)
    Chem.setSMILESWriteKekuleFormParameter(w, kekulized)
    Chem.setSMILESWriteCanonicalFormParameter(w, canonical)
    Chem.setSMILESRecordFormatParameter(w, 'S')
    Chem.setSMILESWriteAtomStereoParameter(w, atom_stereo)
    Chem.setSMILESWriteBondStereoParameter(w, bond_stereo)
    Chem.setSMILESNoOrganicSubsetParameter(w, False)
    Chem.setOrdinaryHydrogenDepleteParameter(w, hydrogen_deplete)
    Chem.calcImplicitHydrogenCounts(mol, True)
    w.write(mol)
    w.close()
    return stream.value


def mol_to_sdf(molecules, path, multiconf=True):
    if not isinstance(molecules, Iterable):
        molecules = [molecules]
    w = Chem.FileSDFMolecularGraphWriter(path)
    Chem.setMultiConfExportParameter(w, multiconf)
    for mol in molecules:
        Chem.calcImplicitHydrogenCounts(mol, False)
        w.write(mol)
    w.close()


def is_macrocyle(mol, ring_size=7):
    """
    Checks if the given molecule contains rings larger than given size. If yes --> macrocyle.
    :param mol:
    :param ring_size:
    :return: Boolean indicating if macrocylce or not
    """
    sssr = Chem.perceiveSSSR(mol)
    if sssr.getSize() > 0:
        if max([frag.getNumAtoms() for frag in sssr]) > ring_size:
            return True
    return False


# def align_molecules_by_pharmacophores(reference: Chem.BasicMolecule=None,
#                                       query: Chem.BasicMolecule=None,
#                                       mode: str="first",
#                                       return_aligned_molecule=False,
#                                       ):
#     """
#     Aligns two molecules with multiple conformations by their pharmacophores.
#     :param reference:
#     :param query:
#     :param mode: one of [first, best]
#     :param return_aligned_molecule:
#     :return:
#     """
#     best_aligment_score = 0
#     best_reference_conformation = None
#     best_tf_matrix = None
#     best_molecule = None
#     for i in range(Chem.getNumConformations(reference)):
#         Chem.applyConformation(reference, i)
#         ph4 = get_pharmacophore(reference)
#         score, tf_matrix, aligned_mol = align_molecule_to_pharmacophore(query=query,
#                                                                         reference=ph4,
#                                                                         mode=mode,
#                                                                         return_aligned_molecule=return_aligned_molecule)
#         if score is not None:
#             if score > best_aligment_score:
#                 best_aligment_score = score
#                 best_tf_matrix = score
#                 best_reference_conformation = i
#                 best_molecule = aligned_mol
#     if not return_aligned_molecule:
#         return best_aligment_score, best_tf_matrix, None
#     else:
#         if best_tf_matrix is not None:
#             Chem.applyConformation(reference, best_reference_conformation)
#             Chem.transform3DCoordinates(reference, best_tf_matrix)
#             return best_aligment_score, best_tf_matrix, best_molecule
#         else:
#             return best_aligment_score, None, None
#
#
# def align_molecule_to_pharmacophore(query: Chem.BasicMolecule=None,
#                                     reference: Pharm.BasicPharmacophore=None,
#                                     mode="first",
#                                     return_aligned_molecule=False,
#                                     ):
#     """
#     Aligns a single molecule with multiple conformations to a pharmacophore.
#     :param query:
#     :param reference:
#     :param mode:
#     :param return_aligned_molecule:
#     :return:
#     """
#
#     if mode not in ["first", "best"]:
#         raise ValueError("Alignment mode not recognized. Needs to be one of [first, best]. "
#                          "%s was given" % mode)
#
#     pharmacophore_aligner = Pharm.PharmacophoreAlignment(True)
#     pharmacophore_aligner.addFeatures(reference, True)
#     scorer = Pharm.PharmacophoreFitScore()
#     best_alignment_score = 0
#     best_tf_matrix = None
#     best_conformation = None
#
#     for i in range(Chem.getNumConformations(query)):
#         Chem.applyConformation(query, i)
#         ph4 = get_pharmacophore(query)
#         pharmacophore_aligner.addFeatures(ph4, False)
#         if mode == "first":
#             if pharmacophore_aligner.nextAlignment():
#                 tf_matrix = pharmacophore_aligner.getTransform()
#                 score = scorer(reference, ph4, tf_matrix)
#                 if score is not None:
#                     if score > best_alignment_score:
#                         best_alignment_score = score
#                         best_tf_matrix = tf_matrix
#                         best_conformation = i
#         elif mode == "best":
#             while pharmacophore_aligner.nextAlignment():
#                 tf_matrix = pharmacophore_aligner.getTransform()
#                 score = scorer(reference, ph4, tf_matrix)
#                 if score is not None:
#                     if score > best_alignment_score:
#                         best_alignment_score = score
#                         best_tf_matrix = tf_matrix
#                         best_conformation = i
#         pharmacophore_aligner.clearEntities(False)
#
#     if not return_aligned_molecule:
#         return best_alignment_score, best_tf_matrix, None
#     else:
#         if best_tf_matrix is not None:
#             Chem.applyConformation(query, best_conformation)
#             Chem.transform3DCoordinates(query, best_tf_matrix)
#             return best_alignment_score, best_tf_matrix, query
#         else:
#             return best_alignment_score, None, None


def calculate_molecule_hashcode(mol, stereo=True):
    Chem.makeHydrogenDeplete(mol)
    Chem.calcImplicitHydrogenCounts(mol, True)
    if stereo:
        Chem.calcAtomStereoDescriptors(mol, True)
        Chem.calcBondStereoDescriptors(mol, True)
        Chem.calcCIPPriorities(mol, True)
        Chem.calcAtomCIPConfigurations(mol, True)
        Chem.calcBondCIPConfigurations(mol, True)
        return Chem.calcHashCode(mol)
    else:
        return Chem.calcHashCode(mol,
                                 atom_flags=Chem.AtomPropertyFlag.TYPE | Chem.AtomPropertyFlag.H_COUNT | Chem.AtomPropertyFlag.FORMAL_CHARGE | Chem.AtomPropertyFlag.AROMATICITY,
                                 bond_flags=Chem.BondPropertyFlag.ORDER | Chem.BondPropertyFlag.TOPOLOGY | Chem.BondPropertyFlag.AROMATICITY
                                 )


def calculateStandardProperties(mol):
    standardProperties = {
        'nrAcceptors': [],
        'nrDonors': [],
        # 'nrRings': [],
        'nrRotBonds': [],
        'molWeight': [],
        'nrHeavyAtoms': [],
        'cLogP': [],
        'TPSA': [],
    }

    try:
        iter(mol)
    except:
        mol = [mol]

    for m in mol:
        Chem.calcTopologicalDistanceMatrix(m, True)

        p = getPharmacophore(m)
        hba, hbd = 0, 0
        for f in p:
            if Pharm.getType(f) == Pharm.FeatureType.H_BOND_ACCEPTOR:
                hba += 1
            elif Pharm.getType(f) == Pharm.FeatureType.H_BOND_DONOR:
                hbd += 1

        standardProperties['nrAcceptors'].append(hba)
        standardProperties['nrDonors'].append(hbd)
        standardProperties['molWeight'].append(Chem.calcExplicitMass(m))
        standardProperties['nrHeavyAtoms'].append(Chem.getHeavyAtomCount(m))
        standardProperties['cLogP'].append(Chem.calcXLogP(m))
        standardProperties['TPSA'].append(Chem.calcTPSA(m))
        standardProperties['nrRotBonds'].append(Chem.getRotatableBondCount(m, False, False))

    return standardProperties
