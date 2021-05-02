import sys
import os
import numpy as np
import CDPL.Chem as Chem
import CDPL.Base as Base
import CDPL.Pharm as Pharm
import CDPL.Math as Math


# RDF parameters
RADIUS = 20
RADINC = 0.1
SEMI = True  # indicate if split or dense RDF codes should be calculated
CLASSIC = True  # indicate if classic RDF should  be calculated, where smoothing factor needs to be calculated
PRIME = False  # indicate if prime number should be used for calculation
FEATURE_TYPES = {
    1: Pharm.FeatureType.HYDROPHOBIC,
    2: Pharm.FeatureType.AROMATIC,
    3: Pharm.FeatureType.NEG_IONIZABLE,
    4: Pharm.FeatureType.POS_IONIZABLE,
    5: Pharm.FeatureType.H_BOND_DONOR,
    6: Pharm.FeatureType.H_BOND_ACCEPTOR,
    7: Pharm.FeatureType.X_VOLUME,
}
FEATURE_TYPES_INVERSE = {value: key for key, value in FEATURE_TYPES.items()}


def save_pharmacophore(pharmacophore: Pharm.BasicPharmacophore, path: str):
    # print("Saving Pharmacophore")
    writer = Pharm.FilePMLFeatureContainerWriter(path)
    writer.write(pharmacophore)
    writer.close()


def load_pml_pharmacophore(path):
    # print("Loading pharmacophore from %s" % path)
    ifs = Base.FileIOStream(path)
    r = Pharm.PMLPharmacophoreReader(ifs)
    pharm = Pharm.BasicPharmacophore()
    try:
        r.read(pharm)
        return pharm
    except:
        return False


def get_interaction_pharmacophore(protein: Chem.BasicMolecule,
                                  ligand: Chem.BasicMolecule,
                                  exclusion_volumes=True,
                                  fuzzy=False,
                                  ) -> Pharm.BasicPharmacophore:
    """
    Create an interaction pharmacophore between a given binding site and a ligand. Both of them need to have matching
    3D coordinates.
    :param protein:
    :param ligand:
    :return:
    """
    assert isinstance(protein, Chem.BasicMolecule) and isinstance(ligand, Chem.BasicMolecule)
    Pharm.prepareForPharmacophoreGeneration(protein)
    Pharm.prepareForPharmacophoreGeneration(ligand)
    Chem.generateHydrogen3DCoordinates(protein, False)
    Chem.generateHydrogen3DCoordinates(ligand, False)

    int_pharm = Pharm.BasicPharmacophore()
    pharm_gen = Pharm.InteractionPharmacophoreGenerator(False, True)  # non-fuzzy core ph4, fuzzy env. ph4
    pharm_gen.addExclusionVolumes(exclusion_volumes)
    pharm_gen.generate(ligand, protein, int_pharm, True)  # True means ligand environment shall be extracted first

    if fuzzy:
        for f in int_pharm:
            if Pharm.getType(f) == 5 or Pharm.getType(f) == 6:
                Pharm.clearOrientation(f)
                Pharm.setGeometry(f, Pharm.FeatureGeometry.SPHERE)
    return int_pharm


def get_pharmacophore(mol: Chem.BasicMolecule, fuzzy=True) -> Pharm.BasicPharmacophore:
    """

    :param mol: Molecule to generate pharmacophore from.
    :param fuzzy: Indicates whether to generate a vector for HBD and HBA or just a sphere. Fuzzy=sphere.
    :return:
    """
    if isinstance(mol, Pharm.BasicPharmacophore):
        return mol
    assert isinstance(mol, Chem.BasicMolecule), "Given object should be of type Chem.BasicMolecule, %s was given" % type(mol)
    Pharm.prepareForPharmacophoreGeneration(mol)  # Fails silently if molecule has coordinates == 0 !!!
    # Chem.makeHydrogenComplete(mol)
    Chem.generateHydrogen3DCoordinates(mol, False)
    pharm = Pharm.BasicPharmacophore()
    pharm_generator = Pharm.DefaultPharmacophoreGenerator(fuzzy)
    pharm_generator.generate(mol, pharm)
    if fuzzy:
        for f in pharm:
            if Pharm.getType(f) == 5 or Pharm.getType(f) == 6:
                Pharm.clearOrientation(f)
                Pharm.setGeometry(f, Pharm.FeatureGeometry.SPHERE)
    return pharm


# def align_multiple_pharmacophores_to_pharmacophore(reference: Pharm.BasicPharmacophore=None,
#                                                    query: list=None,
#                                                    alignment: str ="first",
#                                                    return_transformation_matrix=False,
#                                                    ):
#     """
#     Align the given pharmacophores and retrieve the corresponding alignment scores.
#     :param reference: Pharmacophore to align to
#     :param query: list of pharmacophores which should be aligned to the first
#     :param alignment: 'first' or 'best'. Returns either the first or best, determined by alignment score, pharmacophore
#     found, respectively.
#     :param return_transformation_matrix:
#     """
#     scorer = Pharm.PharmacophoreFitScore()
#     scores = {}
#
#     if not isinstance(query, (list, tuple, set)):
#         query = [query]
#     if reference is None:
#         # align = Pharm.PharmacophoreAlignment(False)
#         if len(query) >= 2:
#             reference = query[0]
#             query = query[1:]
#         else:
#             raise ValueError("List of ph4 to align needs to contain at least two if no reference ph4 is given")
#
#     aligner = Pharm.PharmacophoreAlignment(True)
#     aligner.addFeatures(reference, True)  # add the reference ph4
#
#     for i, ph4 in enumerate(query):
#         if ph4 == reference:
#             scores[ph4] = 1
#             continue
#         # print("Aligning pharmacophore %s of %s" % (str(i), len(to_align)))
#         aligner.addFeatures(ph4, False)
#         if alignment == "first":
#             if aligner.nextAlignment():
#                 tf_matrix = aligner.getTransform()
#                 score = scorer(reference, ph4, tf_matrix)
#                 if return_transformation_matrix:
#                     scores[ph4] = (score, Math.Matrix4D(tf_matrix))  # alignment matrix copy!
#                 else:
#                     scores[ph4] = score
#             else:
#                 scores[ph4] = 0  # no alignment
#         elif alignment == "best":
#             best_score = None
#             best_tf_matrix = None
#             while aligner.nextAlignment():  # skips if alignment is not possible
#                 tf_matrix = aligner.getTransform()
#                 score = scorer(reference, ph4, tf_matrix)
#                 if best_score is not None:
#                     if score > best_score:
#                         best_score = score
#                         best_tf_matrix.assign(tf_matrix)
#                 else:
#                     best_score = score
#                     best_tf_matrix = Math.Matrix4D(tf_matrix)
#             if return_transformation_matrix:
#                 scores[ph4] = (best_score, best_tf_matrix)
#             else:
#                 scores[ph4] = best_score
#         else:
#             raise ValueError("%s not in list of allowed alignment types ['first', 'best']" % alignment)
#         aligner.clearEntities(False)
#
#     if len(scores) == 1:
#         return list(scores.values())[0]
#     return scores
#
#
# def align_pharmacophore_to_pharmacophore(reference,
#                                          query,
#                                          mode="first",  # 'first' or 'best'
#                                          return_transformation_matrix=False):
#     assert isinstance(reference, Pharm.BasicPharmacophore) and isinstance(query, Pharm.BasicPharmacophore), "Given pharmacophores need to be of type: CDPL.Pharm.BasicPharmacophore"
#
#     scorer = Pharm.PharmacophoreFitScore()
#     aligner = Pharm.PharmacophoreAlignment(True)
#     aligner.addFeatures(reference, True)  # add the reference ph4
#     aligner.addFeatures(query, False)
#     if mode == "first":
#         if aligner.nextAlignment():
#             tf_matrix = aligner.getTransform()
#             score = scorer(reference, query, tf_matrix)
#             if return_transformation_matrix:
#                 return score, tf_matrix
#             else:
#                 return score, None
#         else:
#             return None, None
#     elif mode == "best":
#         best_score = None
#         best_matrix = None
#         while aligner.nextAlignment():
#             tf_matrix = aligner.getTransform()
#             score = scorer(reference, query, tf_matrix)
#             if best_score is not None:
#                 if score > best_score:
#                     best_score = score
#                     best_matrix.assign(tf_matrix)
#             else:
#                 best_score = score
#                 best_matrix = Math.Matrix4D(tf_matrix)
#         if return_transformation_matrix:
#             return best_score, best_matrix
#         else:
#             return best_score, None


def calculate_rdf(pharm, radius=20, radInc=0.1, semi=True, classic=True, prime=False):
    """
    Calculate the RDF code for a given pharmacophore.
    :param pharm:
    :param radius:
    :param radInc:
    :param semi: Split or Dense RDF Codes (split if True)
    :param classic: Classic RDFs calculate the smoothing factor
    :param prime: Prime number used for calculation
    :return:
    """
    step = np.long(RADIUS / RADINC)
    entity = False
    if SEMI:
        rdf_calc = Pharm.PharmacophoreRDFDescriptorCalculator()
    else:
        rdf_calc = Pharm.FeatureRDFCodeCalculator()
        entity = True
        rdf_calc.setEntityPairWeightFunction(featurePairWeightFunc)

    smoothing_factor = getSmoothingFactor()
    if not CLASSIC:
        rdf_calc.enableDistanceToIntervalCenterRounding(True)

    rdf_calc.setSmoothingFactor(smoothing_factor)
    rdf_calc.setNumSteps(step - 1)
    rdf_calc.setRadiusIncrement(RADINC)
    if PRIME:
        if entity:
            rdf_calc.setEntityPairWeightFunction(featurePairWeightFunc)
        else:
            rdf_calc.setFeaturePairWeightFunction(featurePairWeightFuncSemi)
    rdf = Math.DVector()
    rdf_calc.calculate(pharm, rdf)
    return np.array(rdf)


def rdf_cos_distance(rdf_1, rdf_2, cosine_sim=True):
    """
    Calculates the cosine ____ between two numpy arrays.
    Cosine distance is defined as 1 - cosine similarity.
    :param rdf_1:
    :param rdf_2:
    :param cosine_sim: If true, calcualtes the cosine similarity instead of cosine distance.
    :return:
    """
    from sklearn.metrics.pairwise import cosine_similarity

    sim = cosine_similarity(rdf_1.reshape(1, -1), rdf_2.reshape(1, -1)).flatten()[0].round(5)
    if cosine_sim:
        return sim
    return 1-sim


def getSmoothingFactor():
    inter = 4 * np.log(2)
    smoothing = np.floor(inter/np.power(float(RADINC)/2, 2))
    return smoothing


def atomPairWeightFuncSemi(atom1, atom2, atm_type):
    at1 = Chem.getType(atom1)
    at2 = Chem.getType(atom2)

    if at1 != atm_type and at2 != atm_type:
        return 0

    if at2 == atm_type:
        return Math.prime(at1)

    return Math.prime(at2)


def featurePairWeightFuncSemi(ftr1, ftr2, ftr_type):
    ft1 = Pharm.getType(ftr1) + 1
    ft2 = Pharm.getType(ftr2) + 1

    if ft1 != ftr_type and ft2 != ftr_type:
        return 0

    if ft2 == ftr_type:
        return Math.prime(ft1)

    return Math.prime(ft2)


def atomPairWeightFunc(atom1, atom2):
    at1 = Chem.getType(atom1)
    at2 = Chem.getType(atom2)

    return at1 * at2
    # if at1 > at2:
    #     return Math.prime((Chem.AtomType.MAX_TYPE + 1) * at2 + at1)

    # return Math.prime((Chem.AtomType.MAX_TYPE + 1) * at1 + at2)


def featurePairWeightFunc(ftr1, ftr2):
    ft1 = Pharm.getType(ftr1) + 1
    ft2 = Pharm.getType(ftr2) + 1

    return ft1 * ft2
    # if ft1 > ft2:
    #     return Math.prime((Pharm.FeatureType.MAX_TYPE + 1) * ft2 + ft1)

    # return Math.prime((Pharm.FeatureType.MAX_TYPE + 1) * ft1 + ft2)


def pharmacophore_coords_to_numpy(ph4, return_feature_type=False):
    coords = []
    if return_feature_type:
        features = []
        for feature in ph4:
            coords.append(Chem.get3DCoordinates(feature))
            features.append(Pharm.getType(feature))
        return coords, features
    else:
        for feature in ph4:
            coords.append(Chem.get3DCoordinates(feature))
        return np.array(coords)


def get_number_of_matching_feature_pairs(reference, query):
    """
    Compare the features of both given pharmacophores and return the number of matching faetures as well as found
    feature pairs.
    :param referece:
    :param query:
    :return:
    """
    tfa = Pharm.TopologicalFeatureAlignment()
    tfa.setEntityMatchFunction(Pharm.FeatureTypeMatchFunctor())
    tfa.setEntityPairMatchFunction(Pharm.FeaturePairDistanceMatchFunctor(False))

    query_mapping = Pharm.FeatureMapping()

    # set reference features
    for feature in reference:
        tfa.addEntity(feature, True)

    # set query features
    for feature in query:
        tfa.addEntity(feature, False)

    possible_alignments = 0
    while tfa.nextAlignment(query_mapping):
        possible_alignments += 1

    return possible_alignments


