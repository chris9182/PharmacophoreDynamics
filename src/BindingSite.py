from abc import ABCMeta, abstractmethod
import os
from typing import Tuple, Any

import numpy as np
import open3d as o3d
import CDPL.Chem as Chem
from scipy import sparse
import scipy.io as sio

from src.utils.Protein import Protein
from src.utils.MathTools import rigid_transform_3D, rotate3DObject, translate3DObject


PointCloud = o3d.geometry.PointCloud
Vector3dVector = o3d.utility.Vector3DVector


# define some aliases
Surface = np.array
Histogram = np.array
Keypoints = np.array
HistogramDistances = np.array
EuclideanDistance = np.array
ContactMatrix = sparse.lil_matrix


class BaseBindingSiteFinder(metaclass=ABCMeta):

    """
    Takes two proteins and returns the two proteins aligned.
    """

    def __init__(self, p1: Protein, p2: Protein, *args, **kwargs):
        self.p1 = p1
        self.p2 = p2

    def __call__(self, nrSolutions: int = 10, *args, **kwargs):
        """
        Find the binding site between
        :param args:
        :param kwargs:
        :return:
        """
        solutions = []
        i = 0
        for solution in self:
            solutions.append(solution)
            i += 1
            if i == nrSolutions:
                return solutions

    def __iter__(self):
        preparedProteins, nextInput = self.prepareProteins(self.p1, self.p2)

        finished = nextInput[0]
        while not finished:
            solution, nextInput = self.findBindingSite(preparedProteins, nextInput)
            finished = nextInput[0]
            yield solution

    @abstractmethod
    def prepareProteins(self, p1, p2, *args, **kwargs) -> Tuple[Tuple[Protein, Protein], Any]:
        """
        Takes two proteins and applies some preparation steps. The output is passed to the findBindinSite function,
        which yields the solutions -> aligned proteins.
        :param p1:
        :param p2:
        :param args:
        :param kwargs:
        :return: Tuple of aligned prepared proteins as well as the next input. The first element of the next input
        should be a boolean indicating whether its finished or not.
        """
        raise NotImplementedError

    @abstractmethod
    def findBindingSite(self, preparedProteins, nextInput, *args, **kwargs) -> Tuple[Tuple[Protein, Protein], Any]:
        """
        Takes two proteins in a random orientation at a random distance to each other as input.
        :param preparedProteins:
        :param nextInput:
        :param args:
        :param kwargs:
        :return:
        """
        raise NotImplementedError


class ShapeMatching(BaseBindingSiteFinder):

    def __init__(self, p1: Protein, p2: Protein, *args, **kwargs):
        self.dataPath = '../data/compute/'
        super(ShapeMatching, self).__init__(p1, p2)

    def prepareProteins(self, p1, p2, *args, **kwargs):
        if Chem.getName(p1) == '':
            Chem.setName(p1, 'p1')
        if Chem.getName(p2) == '':
            Chem.setName(p2, 'p2')

        # calculate molecular surface
        mesh1 = self.calculateMolecularSurface(p1)
        mesh2 = self.calculateMolecularSurface(p2)

        # make point clouds
        pcd1, surface1, keypoints1 = self.makePointCloud(mesh1, counterParty=False)
        pcd2, surface2, keypoints2 = self.makePointCloud(mesh2, counterParty=True)

        # compute histogram features
        pcd1, surface1, histogram1 = self.computePointEncoding(pcd1, surface1, keypoints1)
        pcd2, surface2, histogram2 = self.computePointEncoding(pcd2, surface2, keypoints2)

        # compute distances
        histDistances, eucDist1, eucDist2 = self.computeDistances(histogram1, histogram2, surface1, surface2)

        # compute graph
        pairs, contactMatrix, nrContacts = self.computeGraph(surface1, surface2, histDistances, eucDist1, eucDist2)

        finished = not (nrContacts > 0)
        nextInput = (finished, pairs, contactMatrix, nrContacts, pcd1, pcd2, surface1, surface2)
        return (p1, p2), nextInput

    def findBindingSite(self, preparedProteins, nextInput, *args, **kwargs) -> Tuple[Tuple[Protein, Protein], Any]:
        p1, p2 = preparedProteins
        finished, pairs, contactMatrix, nrContacts, pcd1, pcd2, surface1, surface2 = nextInput

        # computing clique
        clique = self.computeClique(contactMatrix)
        # self.visualizeSolution(clique, pcd1, pcd2, surface1, surface2)

        # set points
        points1, points2 = np.zeros((len(clique), 3)), np.zeros((len(clique), 3))
        index = 0
        for v in clique:
            p = pairs[v]
            points1[index] = surface1[p[0]]
            points2[index] = surface2[p[0]]
            index += 1

        c, R, t = rigid_transform_3D(points1, points2, False)  # get transformation matrices -> kapsch

        alignedP1 = Protein(p1)
        alignedP2 = Protein(p2)
        alignedP1.rotate(R).translate(t)
        alignedP2.rotate(R).translate(t)

        for index in clique:
            for index2 in clique:
                nrContacts -= 1
        finished = not (nrContacts > 0)

        nextInput = (finished, pairs, contactMatrix, nrContacts, pcd1, pcd2, surface1, surface2)
        return (alignedP1, alignedP2), nextInput

    def computeClique(self, contactMatrix):
        sio.mmwrite('{}contact_matrix.mtx'.format(self.dataPath), contactMatrix)
        os.system('OMP_NUM_THREADS=12 ../TEASER-plusplus/build/pmc-build/pmc_main -f {d}/contact_matrix.mtx -a 0 > {d}/result.txt'.format(d=self.dataPath))
        os.system('grep -hr \"Maximum clique\" {d}/result.txt | tail -1 > {d}/clique.txt'.format(d=self.dataPath))
        with open('{}/clique.txt'.format(self.dataPath), 'r') as myfile:
            data = myfile.readlines()
        stringVal = data[0]
        print(stringVal)
        start = stringVal.index(':') + 2
        end = stringVal.index('\n') - 1
        stringVal = stringVal[start:end]
        strArr = stringVal.split()
        intArr = [int(numeric_string) - 1 for numeric_string in strArr]  # minus 1 because scipy outputs with 1 based array
        return intArr

    def calculateMolecularSurface(self,
                                  p: Protein,
                                  h: int = 2,
                                  s: int = 3
                                  ) -> o3d.cpu.pybind.geometry.TriangleMesh:
        fileName = '{}{}'.format(self.dataPath, Chem.getName(p))

        # save proteins as PDB to be read by other programs
        p.toFile(fileName)

        # calculate surface with EDTSurf
        os.system('../EDTSurf/EDTSurf -i {i} -h {h} -s {s}'.format(i='{}.pdb'.format(fileName), h=h, s=s))
        mesh = o3d.io.read_triangle_mesh('{}.ply'.format(fileName))

        # TODO: add cleanup
        return mesh

    def makePointCloud(self, mesh: o3d.cpu.pybind.geometry.TriangleMesh, counterParty=False) -> Tuple[PointCloud, Surface, Keypoints]:
        mesh.compute_vertex_normals()
        pointCloud = PointCloud()
        pointCloud.points = Vector3dVector(mesh.vertices)
        pointCloud.normals = Vector3dVector(
            mesh.vertex_normals if not counterParty else np.negative(mesh.vertex_normals))

        # compute key points
        keypoints = o3d.geometry.keypoint.compute_iss_keypoints(pointCloud)
        keypointsPts = np.array(keypoints.points)

        surface = np.concatenate((np.array(pointCloud.points), keypointsPts), axis=0)
        normals = np.concatenate((np.array(pointCloud.normals), np.array(keypoints.normals)), axis=0)

        finalPointCloud = PointCloud()
        finalPointCloud.points = Vector3dVector(surface)
        finalPointCloud.normals = Vector3dVector(normals)
        return finalPointCloud, surface, keypoints

    def computePointEncoding(self,
                             pointCloud: PointCloud,
                             surface: Surface,
                             keypoints: Keypoints,
                             searchRadius: float = 5,
                             portion: int = 1,
                             ) -> Tuple[PointCloud, Surface, Histogram]:
        features = o3d.pipelines.registration.compute.fpfh_feature(pointCloud, o3d.geometry.KDTreeSearchParamRadius(
            radius=searchRadius))

        numberOfPoints = len(keypoints.points)
        histogram = np.array(features.data).transponse()
        surface = surface[-numberOfPoints:, :]
        # normals = np.array(pointCloud.normals)[-numberOfPoints:, :]
        histogram = histogram[-numberOfPoints:, :]

        choice = np.random.choice(surface.shape[0], int(len(surface) / portion), replace=False)
        surface = surface[choice, :]
        normals = np.array(pointCloud.normals)[choice, :]  # TODO: check shape -> really overwrite from before?
        histogram = histogram[choice, :]

        pointCloud.points = Vector3dVector(surface)
        pointCloud.normals = Vector3dVector(normals)

        # TODO: maybe add pharmacophore encoding here
        return pointCloud, surface, histogram

    def computeDistances(self,
                         histogram1: Histogram,
                         histogram2: Histogram,
                         surface1: Surface,
                         surface2: Surface,
                         ) -> Tuple[HistogramDistances, EuclideanDistance, EuclideanDistance]:
        from sklearn.metrics import euclidean_distances

        histogramDistances = euclidean_distances(histogram1, histogram2)
        eucDist1 = euclidean_distances(surface1)
        eucDist2 = euclidean_distances(surface2)
        return histogramDistances, eucDist1, eucDist2

    def computeGraph(self,
                     surface1: Surface,
                     surface2: Surface,
                     histogramDistances: HistogramDistances,
                     eucDist1: EuclideanDistance,
                     eucDist2: EuclideanDistance,
                     maxDiff: int = 1,
                     maxDist: int = 15,
                     maxHistDiff: int = 30,
                     ) -> Tuple[np.array, ContactMatrix, int]:
        pairs = []
        for i in range(len(surface1)):
            for j in range(len(surface2)):
                histDist = histogramDistances[i, j]
                if histDist < maxHistDiff:
                    pairs.append((i, j))
        pairs2 = np.where(histogramDistances < maxHistDiff)  # TODO debugger -> is it the same? time!

        lengthPairs = len(pairs)
        contactMatrix = sparse.lil_matrix((lengthPairs, lengthPairs), dtype=int)
        nrContacts = 0
        for i in range(lengthPairs):
            p1 = pairs[i]
            for j in range(i+1, lengthPairs):
                p2 = pairs[j]

                if p1[0] == p2[0] or p1[1] == p2[1]:
                    continue

                distance1 = eucDist1[p1[0]][p2[0]]
                if distance1 > maxDist:
                    continue

                distance2 = eucDist2[p1[1]][p2[1]]
                if distance2 > maxDist:
                    continue

                if abs(distance1 - distance2) < maxDiff:
                    contactMatrix[i, j] = 1
                    contactMatrix[j, i] = 1
                    nrContacts += 2
        print('Contact pairs: ', contactMatrix.shape)
        return pairs, contactMatrix, nrContacts

    def visualizeSolution(self,
                          clique: ContactMatrix,
                          pointCloud1: PointCloud,
                          pointCloud2: PointCloud,
                          surface1: Surface,
                          surface2: Surface,
                          ):
        raise NotImplementedError
