import CDPL.Math as Math
import numpy as np


def rotate3DObject(object3D: Math.Vector3DArray, rotMatrix: Math.Matrix3D) -> Math.Vector3DArray:
    object3D.assign(np.matmul(object3D.toArray(False), rotMatrix.toArray()))
    return object3D


def translate3DObject(object3D: Math.Vector3DArray, direction: Math.Vector3D) -> Math.Vector3DArray:
    object3D.assign(object3D.toArray(False) + direction.toArray())
    return object3D


# Implements Kabsch algorithm - best fit.
# Supports scaling (umeyama)
# Compares well to SA results for the same data.
# Input:
#     Nominal  A Nx3 matrix of points
#     Measured B Nx3 matrix of points
# Returns s,R,t
# s = scale B to A
# R = 3x3 rotation matrix (B to A)
# t = 3x1 translation vector (B to A)
# from https://gist.github.com/oshea00/dfb7d657feca009bf4d095d4cb8ea4be
def rigid_transform_3D(A, B, scale):
    assert len(A) == len(B)

    N = A.shape[0]  # total points

    centroid_A = np.mean(A, axis=0)
    centroid_B = np.mean(B, axis=0)

    # center the points
    AA = A - np.tile(centroid_A, (N, 1))
    BB = B - np.tile(centroid_B, (N, 1))

    # dot is matrix multiplication for array
    if scale:
        H = np.transpose(BB) @ AA / N
    else:
        H = np.transpose(BB) @ AA

    U, S, Vt = np.linalg.svd(H)

    R = Vt.T @ U.T

    # special reflection case
    if np.linalg.det(R) < 0:
        print("Reflection detected")
        Vt[2, :] *= -1
        R = Vt.T @ U.T

    if scale:
        varA = np.var(A, axis=0).sum()
        c = 1 / (1 / varA * np.sum(S))  # scale factor
        t = -R @ (centroid_B.T * c) + centroid_A.T
    else:
        c = 1
        t = -R @ centroid_B.T + centroid_A.T

    return c, R, t