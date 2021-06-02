import math

import CDPL.Base as Base
import CDPL.Chem as Chem
import CDPL.Biomol as Biomol
import CDPL.Math as Math
import CDPL.Pharm as Pharm
import pyvdwsurface as surface

import open3d as o3d
import numpy as np
import sklearn
from matplotlib import pyplot as plt
from sklearn.cluster import DBSCAN
from sklearn.cluster import KMeans
from sklearn import preprocessing
from sklearn.decomposition import PCA
from sklearn.metrics import euclidean_distances
from sklearn.metrics.pairwise import cosine_distances
from sklearn_extra.cluster import KMedoids

import MathTools
import ProteinTools
import PharmacophoreTools as PhaTools
from Protein import Protein
from simple_xmeans import XMeans
import networkx as nx
from networkx.algorithms.approximation import clique
from queue import PriorityQueue
from scipy import sparse
import scipy.io as sio
import pcl

import os
import sys


def cosine_similarity_np(list_1, list_2):
    cos_sim = np.dot(list_1, list_2) / (np.norm(list_1) * np.norm(list_2))
    return cos_sim


def custom_draw_geometry(pcdarr):
    vis = o3d.visualization.Visualizer()
    vis.create_window()
    vis.get_render_option().point_size = 15
    for pcd in pcdarr:
        vis.add_geometry(pcd)

    vis.run()
    vis.destroy_window()


def create_spheres(data, color, radius=0.5):
    """
    Create a list of spheres from a 2D numpy array

    Numpy array needs to be N-by-3
    """
    vis_list = []
    for row in data:
        c_pt = row
        mesh_sphere = o3d.geometry.TriangleMesh.create_sphere(radius=radius)
        mesh_sphere.compute_vertex_normals()
        mesh_sphere.paint_uniform_color(color)
        mesh_sphere.translate(c_pt)
        vis_list.append(mesh_sphere)
    return vis_list


def create_pha_spheres(pha, radius=0.5):
    colors = [[1, 1, 1], [1, 1, 1], [1, 1, 1], [0, 1, 1], [0, 1, 0], [1, 0, 1], [1, 1, 0], [1, 1, 1], [1, 1, 1],
              [1, 1, 1], [1, 1, 1], [1, 1, 1]]
    vis_list = []
    for feature in pha:
        featureType = Pharm.getType(feature)
        featureCoords = np.array(Chem.get3DCoordinates(feature))
        mesh_sphere = o3d.geometry.TriangleMesh.create_sphere(radius=radius)
        mesh_sphere.compute_vertex_normals()
        mesh_sphere.paint_uniform_color(np.array(colors[featureType]))
        mesh_sphere.translate(featureCoords)
        vis_list.append(mesh_sphere)
    return vis_list


def encodePhaInfo(surface, pha, invert=False):
    types = [3, 4, 5, 6]
    invertedTypes = [4, 3, 6, 5]
    encoding = np.zeros((len(surface), len(types)))
    for feature in pha:
        featureType = Pharm.getType(feature)
        if featureType not in types:
            continue
        featureCoords = np.array(Chem.get3DCoordinates(feature))
        for i in range(len(surface)):
            pt = surface[i]
            dist = np.linalg.norm(pt - featureCoords)
            if invert:
                index = invertedTypes.index(featureType)
            else:
                index = types.index(featureType)
            encoding[i][index] = max(encoding[i][index], 1 / (1 + dist))
    return encoding


def encodePhaInfo2(surface, pha, invert=False):
    types = [-1, -1, -1, 0, 1, 2, 3, -1, -1, -1, -1, -1]
    invertedTypes = [-1, -1, -1, 1, 0, 3, 2, -1, -1, -1, -1, -1]
    typeCount = 4
    encoding = np.full((len(surface), typeCount), np.inf)
    count = 0
    for feature in pha:
        count = count + 1
        featureType = Pharm.getType(feature)
        if invert:
            index = invertedTypes[featureType]
        else:
            index = types[featureType]
        if index < 0:
            continue
        featureCoords = np.array(Chem.get3DCoordinates(feature))
        for i in range(len(surface)):
            pt = surface[i]
            dist = np.linalg.norm(pt - featureCoords)
            encoding[i][index] = min(encoding[i][index], dist)
    print(count)
    for enc in encoding:
        minV = 0
        for i in range(typeCount):
            if enc[minV] > enc[i]:
                minV = i
        # minDist = enc[minV]
        for i in range(typeCount):
            enc[i] = 0
        # if minDist < 20:
        enc[minV] = 1
    return encoding


if __name__ == '__main__':
    # tf.config.list_physical_devices('GPU')
    # exit()
    print("loading molecules")

    # name1 = "1K74_l_b"
    # name2 = "1K74_r_b"
    # name1 = "1KTZ_l_b"
    # name2 = "1KTZ_r_b"
    name1 = "1MAH_l_b"
    name2 = "1MAH_r_b"

    mol1 = Protein()
    mol1.fromFile("structures/" + name1 + ".pdb")
    mol2 = Protein()
    mol2.fromFile("structures/" + name2 + ".pdb")
    mol1.prepareProteins()
    mol2.prepareProteins()
    Chem.makeHydrogenDeplete(mol1)
    Chem.makeHydrogenDeplete(mol2)
    Biomol.FilePDBMolecularGraphWriter("compute/" + name1 + ".pdb").write(mol1)
    Biomol.FilePDBMolecularGraphWriter("compute/" + name2 + ".pdb").write(mol2)
    # ProteinTools.writePDB("compute/" + name1 + ".pdb", mol1)
    # ProteinTools.writePDB("compute/" + name2 + ".pdb", mol2)

    print("computing molecule surface")

    os.system("~/EDTSurf/EDTSurf -i " + "compute/" + name1 + ".pdb" + " -h 2 -s 3")
    os.system("~/EDTSurf/EDTSurf -i " + "compute/" + name2 + ".pdb" + " -h 2 -s 3 ")  # -p 2.0

    '''
    atomCoords1 = Math.Vector3DArray()
    atomCoords2 = Math.Vector3DArray()

    Chem.get3DCoordinates(mol1, atomCoords1)
    Chem.get3DCoordinates(mol2, atomCoords2)

    
    symbols1 = []
    for atom in mol1.atoms:
        symbols1.append(Chem.getSymbolForType(atom))
    symbols2 = []
    for atom in mol2.atoms:
        symbols2.append(Chem.getSymbolForType(atom))

    npCoords1 = np.array(atomCoords1, dtype=float)
    npCoords2 = np.array(atomCoords2, dtype=float)

    print("computing surface")
    # non deterministic!
    surface1 = surface.vdwsurface(npCoords1, symbols1, density=0.2)
    surface2 = surface.vdwsurface(npCoords2, symbols2, density=0.2)
    surface1np = np.array(surface1)
    surface2np = np.array(surface2)

    pcd1 = o3d.geometry.PointCloud()
    pcd1.points = o3d.utility.Vector3dVector(surface1np)

    pcd2 = o3d.geometry.PointCloud()
    pcd2.points = o3d.utility.Vector3dVector(surface2np)

    print(surface1np.shape)
    print(surface2np.shape)

    pcd1.paint_uniform_color([0.0, 0.0, 1.0])  # show pcd1 in blue
    pcd2.paint_uniform_color([1.0, 0.0, 0.0])  # show pcd2 in red
    # o3d.visualization.draw_geometries([pcd1, pcd2])
    # exit()
    '''
    mesh1 = o3d.io.read_triangle_mesh("compute/" + name1 + ".ply")
    mesh2 = o3d.io.read_triangle_mesh("compute/" + name2 + ".ply")
    mesh1.compute_vertex_normals()
    mesh2.compute_vertex_normals()
    pcd1 = o3d.geometry.PointCloud()
    pcd2 = o3d.geometry.PointCloud()
    pcd1.points = o3d.utility.Vector3dVector(mesh1.vertices)
    pcd1.normals = o3d.utility.Vector3dVector(mesh1.vertex_normals)
    pcd2.points = o3d.utility.Vector3dVector(mesh2.vertices)
    pcd2.normals = o3d.utility.Vector3dVector(np.negative(mesh2.vertex_normals))
    pcd1.paint_uniform_color([0.0, 0.0, 1.0])  # show pcd1 in blue
    pcd2.paint_uniform_color([1.0, 0.0, 0.0])  # show pcd2 in red
    print("estimate normals")
    # pcd1.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamRadius(radius=1))
    # pcd2.estimate_normals(search_param=o3d.geometry.KDTreeSearchParamRadius(radius=1))
    # normalsPlaneParam = 10
    # pcd1.orient_normals_consistent_tangent_plane(normalsPlaneParam)
    # pcd2.orient_normals_consistent_tangent_plane(normalsPlaneParam)
    # pcd1.normalize_normals()
    # pcd2.normalize_normals()
    colors1 = [[0, 0, 1] for i in range(len(pcd1.points))]
    colors2 = [[1, 0, 0] for i in range(len(pcd2.points))]
    pcd1.colors = o3d.utility.Vector3dVector(colors1)
    pcd2.colors = o3d.utility.Vector3dVector(colors2)
    # o3d.visualization.draw_geometries([pcd1, pcd2])
    # exit()

    kps1 = o3d.geometry.keypoint.compute_iss_keypoints(pcd1)
    kps2 = o3d.geometry.keypoint.compute_iss_keypoints(pcd2)
    kps1pts = np.array(kps1.points)
    kps2pts = np.array(kps2.points)

    surface1np = np.concatenate((np.array(pcd1.points), kps1pts), axis=0)
    surface2np = np.concatenate((np.array(pcd2.points), kps2pts), axis=0)

    normals1 = np.concatenate((np.array(pcd1.normals), np.array(kps1.normals)), axis=0)
    normals2 = np.concatenate((np.array(pcd2.normals), np.array(kps2.normals)), axis=0)
    pcd1s = o3d.geometry.PointCloud()
    pcd2s = o3d.geometry.PointCloud()
    pcd1s.points = o3d.utility.Vector3dVector(surface1np)
    pcd2s.points = o3d.utility.Vector3dVector(surface2np)
    pcd1s.normals = o3d.utility.Vector3dVector(normals1)
    pcd2s.normals = o3d.utility.Vector3dVector(normals2)

    # pcd1 = o3d.geometry.keypoint.compute_iss_keypoints(pcd1)
    # pcd2 = o3d.geometry.keypoint.compute_iss_keypoints(pcd2)

    print("computing point encoding")
    # rad = 3
    rad = 2
    '''
    radLocal = 0.5
    chainLength = int(rad / radLocal)
    f1 = o3d.pipelines.registration.compute_fpfh_feature(
        pcd1s, o3d.geometry.KDTreeSearchParamNNChain(radiusLocal=radLocal, chainLength=chainLength))
    f2 = o3d.pipelines.registration.compute_fpfh_feature(
        pcd2s, o3d.geometry.KDTreeSearchParamNNChain(radiusLocal=radLocal, chainLength=chainLength))
    '''
    f1 = o3d.pipelines.registration.compute_fpfh_feature(pcd1s, o3d.geometry.KDTreeSearchParamRadius(radius=rad))
    f2 = o3d.pipelines.registration.compute_fpfh_feature(pcd2s, o3d.geometry.KDTreeSearchParamRadius(radius=rad))

    hist1 = np.array(f1.data).transpose()
    print(hist1.shape)
    hist2 = np.array(f2.data).transpose()
    print(hist2.shape)

    surface1np = surface1np[-len(kps1.points):, :]
    surface2np = surface2np[-len(kps2.points):, :]
    normals1 = np.array(pcd1.normals)[-len(kps1.points):, :]
    normals2 = np.array(pcd2.normals)[-len(kps2.points):, :]
    hist1 = hist1[-len(kps1.points):, :]
    hist2 = hist2[-len(kps2.points):, :]

    portion = 1
    choice1 = np.random.choice(surface1np.shape[0], int(len(surface1np) / portion), replace=False)
    choice2 = np.random.choice(surface2np.shape[0], int(len(surface2np) / portion), replace=False)
    surface1np = surface1np[choice1, :]
    surface2np = surface2np[choice2, :]
    normals1 = np.array(pcd1.normals)[choice1, :]
    normals2 = np.array(pcd2.normals)[choice2, :]
    hist1 = hist1[choice1, :]
    hist2 = hist2[choice2, :]

    pcd1s.points = o3d.utility.Vector3dVector(surface1np)
    pcd2s.points = o3d.utility.Vector3dVector(surface2np)
    pcd1s.normals = o3d.utility.Vector3dVector(normals1)
    pcd2s.normals = o3d.utility.Vector3dVector(normals2)
    print(surface1np.shape)
    print(surface2np.shape)

    # surface1np = np.array(pcd1.points)
    # surface2np = np.array(pcd2.points)

    '''
    print("computing pha enc")
    pha1 = PhaTools.getPharmacophore(mol1)
    pha2 = PhaTools.getPharmacophore(mol2)
    surfacePhaEncoding1 = encodePhaInfo2(surface1np, pha1)
    surfacePhaEncoding2 = encodePhaInfo2(surface2np, pha2, True)
    spheresp1 = create_pha_spheres(pha1)
    spheresp2 = create_pha_spheres(pha2)
    '''

    print("computing distances")
    distancesHist = euclidean_distances(hist1, hist2)
    eucDist1 = euclidean_distances(surface1np)
    eucDist2 = euclidean_distances(surface2np)

    print(hist1[2])
    print(hist2[100])
    print(distancesHist[2][100])
    maxDiff = 1
    maxDist = 15
    maxHistDiff = 30

    # maxDiff = 2
    # maxDist = 15
    # maxHistDiff = 65

    print("computing graph")
    distleni = len(surface1np)
    distlenj = len(surface2np)
    pairs = []
    for i in range(distleni):
        # enc1 = surfacePhaEncoding1[i]
        for j in range(distlenj):
            histDist = distancesHist[i][j]
            # enc2 = surfacePhaEncoding2[j]
            if histDist < maxHistDiff:  # np.array_equal(enc1, enc2):
                pairs.append((i, j))

    distlenp = len(pairs)
    mat = sparse.lil_matrix((distlenp, distlenp), dtype=int)
    for i in range(distlenp):
        p1 = pairs[i]
        for j in range(i + 1, distlenp):
            p2 = pairs[j]
            if p1[0] == p2[0] or p1[1] == p2[1]:
                continue
            dist1 = eucDist1[p1[0]][p2[0]]
            if dist1 > maxDist:
                continue
            dist2 = eucDist2[p1[1]][p2[1]]
            if dist2 > maxDist:
                continue
            if abs(dist1 - dist2) < maxDiff:
                mat[i, j] = 1
                mat[j, i] = 1
        b = ("computing row: " + str(i + 1) + " of " + str(distlenp))
        sys.stdout.write('\r' + b)

    print(mat.shape)
    fillIn = np.count_nonzero(mat.toarray())

    rad = 2
    while fillIn > 0:
        print(str(fillIn) + " out of " + str(distlenp * distlenp) + " = " + str(fillIn / (distlenp * distlenp)))
        print("writing matrix")
        sio.mmwrite("compute/sparse_matrix.mtx", mat)
        print("computing clique")
        os.system("OMP_NUM_THREADS=12 ~/git/pmc/build/pmc_main -f compute/sparse_matrix.mtx -a 0 > compute/result.txt")
        os.system("grep -hr \"Maximum clique\" compute/result.txt | tail -1 > compute/clique.txt")
        with open("compute/clique.txt", "r") as myfile:
            data = myfile.readlines()
        stringVal = data[0]
        print(stringVal)
        start = stringVal.index(':') + 2
        end = stringVal.index('\n') - 1
        stringVal = stringVal[start:end]
        strArr = stringVal.split()
        intArr = [int(numeric_string) - 1 for numeric_string in
                  strArr]  # minus 1 because scipy outputs with 1 based array
        # rows, cols = np.where(mat == 1)
        # edges = zip(rows.tolist(), cols.tolist())
        # gr = nx.Graph()
        # gr.add_edges_from(edges)
        # clique = clique.max_clique(gr)
        print(len(intArr))

        print("setting colors")
        sphere1 = create_spheres(np.array(pcd1s.points), [0, 0, 1], rad / 2)
        sphere2 = create_spheres(np.array(pcd2s.points), [1, 0, 0], rad / 2)
        sphere1p = []
        sphere2p = []
        points1 = np.empty((len(intArr), 3))
        points2 = np.empty((len(intArr), 3))
        ind = 0
        for v in intArr:
            p = pairs[v]  # minus 1 because scipy outputs with 1 based array
            points1[ind] = surface1np[p[0]]
            points2[ind] = surface2np[p[1]]
            sphere1[p[0]].paint_uniform_color([0, 0, 0])
            sphere2[p[1]].paint_uniform_color([0, 0, 0])
            sphere1p.append(sphere1[p[0]])
            sphere2p.append(sphere2[p[1]])
            ind = ind + 1

        visList = [pcd1, pcd2]
        visList.extend(sphere1p)
        visList.extend(sphere2p)
        # visList.extend(spheresp1)
        # visList.extend(spheresp2)
        o3d.visualization.draw_geometries(visList)
        # custom_draw_geometry([pcd1, pcd2,pcd1s,pcd2s])

        c, R, t = MathTools.rigid_transform_3D(points1, points2, False)

        pcd2.rotate(R, center=(0, 0, 0))
        pcd2s.rotate(R, center=(0, 0, 0))

        pcd2.translate(t)
        pcd2s.translate(t)
        surface2np = np.array(pcd2s.points)

        sphere1 = create_spheres(np.array(pcd1s.points), [0, 0, 1], rad / 2)
        sphere2 = create_spheres(np.array(pcd2s.points), [1, 0, 0], rad / 2)
        sphere1p = []
        sphere2p = []
        ind = 0
        for v in intArr:
            p = pairs[v]  # minus 1 because scipy outputs with 1 based array
            sphere1[p[0]].paint_uniform_color([0, 0, 0])
            sphere2[p[1]].paint_uniform_color([0, 0, 0])
            sphere1p.append(sphere1[p[0]])
            sphere2p.append(sphere2[p[1]])
            ind = ind + 1

        visList = [pcd1, pcd2]
        visList.extend(sphere1p)
        visList.extend(sphere2p)
        o3d.visualization.draw_geometries(visList)

        for ind in intArr:
            for ind2 in intArr:
                mat[ind, ind2] = 0
        fillIn = np.count_nonzero(mat.toarray())
    '''
    pcd1 = o3d.geometry.PointCloud()
    pcd2 = o3d.geometry.PointCloud()
    pcd1.points = o3d.utility.Vector3dVector(surface1np)
    pcd2.points = o3d.utility.Vector3dVector(surface2np)

    pcd1 = o3d.geometry.PointCloud()
    pcd2 = o3d.geometry.PointCloud()
    pcd1.points = o3d.utility.Vector3dVector(hist1[:, :3])
    pcd2.points = o3d.utility.Vector3dVector(hist2[:, :3])
    
    o3d.visualization.draw_geometries([pcd1, pcd2])
    '''
    '''
    print("computing pharmacophore encoding")

    pha1 = PhaTools.getPharmacophore(mol1)
    pha2 = PhaTools.getPharmacophore(mol2)
    surfacePhaEncoding1 = encodePhaInfo(surface1np, pha1)
    surfacePhaEncoding2 = encodePhaInfo(surface2np, pha2, True)

    data1 = np.concatenate(((hist1, surfacePhaEncoding1)), axis=1)
    data2 = np.concatenate(((hist2, surfacePhaEncoding2)), axis=1)
    print(data1.shape)
    print(data2.shape)
    '''

    '''
    print(hist1[0])

    print("computing point correspondence")
    distances = cosine_distances(hist1, hist2)
    distleni = len(distances)
    distlenj = len(distances[0])
    q = PriorityQueue()
    for i in range(distleni):
        for j in range(distlenj):
            dist = distances[i][j]
            q.put((dist, i, j))
    count = 0
    while True:
        count = count + 1
        print(count)
        min = q.get()
        first = min[1]
        second = min[2]
        print(min[0])

        colors1 = [[0, 0, 1] for i in range(len(pcd1.points))]
        colors2 = [[1, 0, 0] for i in range(len(pcd2.points))]
        colors1[first] = [0, 0, 0]
        colors2[second] = [0, 0, 0]
        pcd1.colors = o3d.utility.Vector3dVector(colors1)
        pcd2.colors = o3d.utility.Vector3dVector(colors2)
        # if np.linalg.norm(surface1np[first] - surface2np[second]) < 10:
        custom_draw_geometry([pcd1, pcd2])
        distances[first][second] = 1000
    '''
    '''
    transformed1 = hist1
    transformed2 = hist2
    pointData = np.append(transformed1, transformed2, axis=0)

    
    print("computing pharmacophore encoding")

    pha1 = PhaTools.getPharmacophore(mol1)
    pha2 = PhaTools.getPharmacophore(mol2)
    surfacePhaEncoding1 = encodePhaInfo(surface1np, pha1)
    surfacePhaEncoding2 = encodePhaInfo(surface2np, pha2, True)
    phaData = np.append(surfacePhaEncoding1, surfacePhaEncoding2, axis=0)
    # phaData = preprocessing.normalize(phaData, norm='max')
    # print(phaData[0])

    print("collecting and transforming")
    data1 = np.concatenate(((transformed1, surfacePhaEncoding1)), axis=1)
    data2 = np.concatenate(((transformed2, surfacePhaEncoding2)), axis=1)
    clusterData = np.append(data1, data2, axis=0)
    # pca = PCA(n_components=10)
    # clusterData = pca.fit_transform(clusterData)
    print(clusterData.shape)
    

    print("clustering")
    db = DBSCAN(eps=0.05, min_samples=10, metric='cosine').fit(pointData)
    # db = DBSCAN(eps=0.1, min_samples=10, metric='euclidean').fit(pointData)
    # core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    # core_samples_mask[db.core_sample_indices_] = True

    #db = KMedoids(n_clusters=20, random_state=0, metric='cosine').fit(pointData)
    labels = db.labels_
    max_label = labels.max()
    print(f"point cloud has {max_label + 1} clusters")
    colors = plt.get_cmap("tab20")(labels / (max_label if max_label > 0 else 1))
    colors[labels < 0] = 0

    borderIndex = len(transformed1)

    pcd1.colors = o3d.utility.Vector3dVector(colors[:borderIndex, :3])
    pcd2.colors = o3d.utility.Vector3dVector(colors[borderIndex:, :3])
    # pcd1.paint_uniform_color([0.0, 0.0, 1.0])  # show pcd1 in blue
    # pcd2.paint_uniform_color([1.0, 0.0, 0.0])  # show pcd2 in red

    o3d.visualization.draw_geometries([pcd1, pcd2])

    pcd1.points = o3d.utility.Vector3dVector(pointData[:borderIndex, :3])
    pcd2.points = o3d.utility.Vector3dVector(pointData[borderIndex:, :3])

    o3d.visualization.draw_geometries([pcd1, pcd2])
    '''
