import CDPL.Math as Math
import numpy as np


def rotate3DObject(object3D: Math.Vector3DArray, rotMatrix: Math.Matrix3D) -> Math.Vector3DArray:
    object3D.assign(np.matmul(object3D.toArray(False), rotMatrix.toArray()))
    return object3D


def translate3DObject(object3D: Math.Vector3DArray, direction: Math.Vector3D) -> Math.Vector3DArray:
    object3D.assign(object3D.toArray(False) + direction.toArray())
    return object3D
