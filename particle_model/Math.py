""" Module containing pure mathematic operations """

from particle_model.Debug import profile, logger
import scipy.linalg as la
import numpy

def invert(mat):
    """ Hard coded 2D matrix inverse."""
    if mat.shape == (2, 2):
        return (numpy.array(((mat[1, 1], -mat[0, 1]),
                             (-mat[1, 0], mat[0, 0])))
                /(mat[0, 0]*mat[1, 1]-mat[0, 1]*mat[1, 0]))

    # otherwise use numpy
    return la.inv(mat)

def cross(vec1, vec2):
    """Return cross product of 3-tuples x and y."""
    out = numpy.zeros(3)
    out[0] = vec1[1]*vec2[2]-vec1[2]*vec2[1]
    out[1] = vec1[2]*vec2[0]-vec1[0]*vec2[2]
    out[2] = vec1[0]*vec2[1]-vec1[1]*vec2[0]
    return out
