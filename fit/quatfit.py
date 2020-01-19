# coding=utf-8
from __future__ import print_function

from math import fabs, sqrt

from misc import transpose


##############################################################################
# The program to superimpose atoms of two molecules by quaternion method
#
# David J. Heisterberg
# The Ohio Supercomputer Center
# 1224 Kinnear Rd.
# Columbus, OH  43212-1163
# (614)292-6036
# djh@ccl.net    djh@ohstpy.bitnet    ohstpy::djh
#
# Translated to C from fitest.f program and interfaced with Xmol program
# by Jan Labanowski,  jkl@ccl.net   jkl@ohstpy.bitnet   ohstpy::jkl
#
# Translated to python from quatfit.c
# by Thomas J. L. Mustard, mustardt@onid.orst.edu
#
# Copyright: Ohio Supercomputer Center, David J. Heisterberg, 1990.
# The program can be copied and distributed freely, provided that
# this copyright in not removed. You may acknowledge the use of the
# program in published material as:
# David J. Heisterberg, 1990, unpublished results.
#


class Atom(object):
    def __init__(self):
        # Hold the string for the element
        self.e = ''
        # Hold the int for the element
        self.en = int(0)
        # Hold the float for the x coordinate
        self.x = float(0)
        # Hold the float for the y coordinate
        self.y = float(0)
        # Hold the float for the z coordinate
        self.z = float(0)
        # Hold the float for the x coordinate
        self.mx = ''
        # Hold the float for the y coordinate
        self.my = ''
        # Hold the float for the z coordinate
        self.mz = ''


def createRefGeom(reflol, fitlol, pairs):
    """
    ### --- Create the reference XYZ files for future alignment --- ###
    :param reflol:
    :param fitlol:
    :param pairs:
    :return:
    """
    # Create ref_xyz and fit_xyz
    ref_xyz = []  # [0] * len(pairs)
    fit_xyz = []  # [0] * len(pairs)
    # Copy the pair atoms from the original molecule into the reference xyz file
    for i, _ in enumerate(pairs):
        # Copy the ref pairs
        ref_xyz.append(Atom())
        ref_xyz[i].e = reflol[pairs[i][0] - 1].e
        ref_xyz[i].x = reflol[pairs[i][0] - 1].x
        ref_xyz[i].y = reflol[pairs[i][0] - 1].y
        ref_xyz[i].z = reflol[pairs[i][0] - 1].z
        ref_xyz[i].mx = reflol[pairs[i][0] - 1].mx
        ref_xyz[i].my = reflol[pairs[i][0] - 1].my
        ref_xyz[i].mz = reflol[pairs[i][0] - 1].mz
        # Copy the fit pairs
        fit_xyz[i] = Atom()
        fit_xyz[i].e = fitlol[pairs[i][1] - 1].e
        fit_xyz[i].x = fitlol[pairs[i][1] - 1].x
        fit_xyz[i].y = fitlol[pairs[i][1] - 1].y
        fit_xyz[i].z = fitlol[pairs[i][1] - 1].z
        fit_xyz[i].mx = fitlol[pairs[i][1] - 1].mx
        fit_xyz[i].my = fitlol[pairs[i][1] - 1].my
        fit_xyz[i].mz = fitlol[pairs[i][1] - 1].mz
    return ref_xyz, fit_xyz


def print_atom(ofilelol):
    """
    ### --- Print the molecule class object --- ###
    :param ofilelol:
    :return:
    """
    for i in range(len(ofilelol)):
        print(ofilelol[i].e + "  " + str("{:.6f}".format(ofilelol[i].x)) + "  " + str(
            "{:.6f}".format(ofilelol[i].y)) + "  " + str("{:.6f}".format(ofilelol[i].z)) + "  " + str(
            "{:.6f}".format(ofilelol[i].charge)) + "  " + str("{:.6f}".format(ofilelol[i].mx)) + "  " + str(
            "{:.6f}".format(ofilelol[i].my)) + "  " + str("{:.6f}".format(ofilelol[i].mz)))
    return


def center(filelol, centerswitch, centerxyz):
    """
    ### --- center the coordinates, or translate them to some xyz --- ###
    CENTER
     center or translate a molecule.
     atomnum (n) - number of atoms
     filelol (x) - on input  - original xyz coordinates of a molecule
         on output - moved xyz coordinates (see io for modes).

     centerswitch (io) - 1 weighted geometric center of the molecule will be at (0,0,0)
          2 molecule will be moved by a vector -center (i.e., components of a vector center
            will be subtracted from atom coordinates).
          3 molecule will be moved by a vector +center (i.e., components of a vector center
            will be added atom coordinates).

     centerxyz (o) - if centerswitch=1, output, center of original coordinates
         if centerswitch=2, input, vector center will be subtracted from atomic coordinates
         if centerswitch=3, input, vector center will be added to atomic coordinates

    """
    modif = float(0.00)
    # int i
    if centerswitch == 2:
        modif = -1.0
    elif centerswitch == 3:
        modif = 1.0
    else:
        modif = -1.0
        centerxyz[0] = float(0.0)
        centerxyz[1] = float(0.0)
        centerxyz[2] = float(0.0)
        for i in range(len(filelol)):
            centerxyz[0] += filelol[i][0]
            centerxyz[1] += filelol[i][1]
            centerxyz[2] += filelol[i][2]
    for i in range(len(filelol)):
        filelol[i][0] += modif * centerxyz[0]
        filelol[i][1] += modif * centerxyz[1]
        filelol[i][2] += modif * centerxyz[2]
    return centerxyz, filelol


def rotmol(filelol, rotmat):
    """
    ROTMOL
    rotate a molecule
    n - number of atoms
    filelol (x) - input coordinates
    filelol (y) - rotated coordinates y = u * x
    rotmat (u) - left rotation matrix
    """
    yx = float(0.0)
    yy = float(0.0)
    yz = float(0.0)
    for i in range(len(filelol)):
        yx = rotmat[0][0] * filelol[i][0] + rotmat[1][0] * filelol[i][1] + rotmat[2][0] * filelol[i][2]
        yy = rotmat[0][1] * filelol[i][0] + rotmat[1][1] * filelol[i][1] + rotmat[2][1] * filelol[i][2]
        yz = rotmat[0][2] * filelol[i][0] + rotmat[1][2] * filelol[i][1] + rotmat[2][2] * filelol[i][2]
        filelol[i][0] = yx  # x
        filelol[i][1] = yy  # y
        filelol[i][2] = yz  # z
    return filelol


def jacobi(matrix, maxsweeps):
    """
    JACOBI
    Jacobi diagonalizer with sorted output. It is only good for 4x4 matrices.
    (was too lazy to do pointers...)
    matrix (a) - input: matrix to diagonalize
    eigenvect (v) - output: eigenvectors
    eigenval (d) - output: eigenvalues
    maxsweeps (nrot) - input: maximum number of sweeps
    """
    eigenvect = [[float(0.0) for x in range(4)] for x in range(4)]
    eigenval = [float(0.0) for x in range(4)]
    onorm = float(0.0)
    dnorm = float(0.0)
    b = float(0.0)
    dma = float(0.0)
    q = float(0.0)
    t = float(0.0)
    c = float(0.0)
    s = float(0.0)
    atemp = float(0.0)
    vtemp = float(0.0)
    dtemp = float(0.0)
    m = 0
    for j in range(4):
        # for i in range(4):
        #  eigenvect[i][j] = 0.0
        eigenvect[j][j] = 1.0
        eigenval[j] = matrix[j][j]
    for m in range(maxsweeps):
        dnorm = 0.0
        onorm = 0.0
        for j in range(4):
            dnorm += fabs(eigenval[j])
            for i in range(j):
                onorm += fabs(matrix[i][j])
        if onorm / dnorm <= 1.0e-12: break  # goto Exit_now;
        for j in range(1, 4):
            for i in range(j):
                b = matrix[i][j]
                if fabs(b) > 0.0:
                    dma = eigenval[j] - eigenval[i]
                    if (fabs(dma) + fabs(b)) <= fabs(dma):
                        t = b / dma
                    else:
                        q = 0.5 * dma / b
                        t = 1.0 / (fabs(q) + sqrt(1.0 + q * q))
                        if q < 0.0:
                            t = -t
                    c = 1.0 / sqrt(t * t + 1.0)
                    s = t * c
                    matrix[i][j] = 0.0
                    for k in range(i):
                        atemp = c * matrix[k][i] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[k][i] + c * matrix[k][j]
                        matrix[k][i] = atemp
                    for k in range(i + 1, j):
                        atemp = c * matrix[i][k] - s * matrix[k][j]
                        matrix[k][j] = s * matrix[i][k] + c * matrix[k][j]
                        matrix[i][k] = atemp
                    for k in range(j + 1, 4):
                        atemp = c * matrix[i][k] - s * matrix[j][k]
                        matrix[j][k] = s * matrix[i][k] + c * matrix[j][k]
                        matrix[i][k] = atemp
                    for k in range(4):
                        vtemp = c * eigenvect[k][i] - s * eigenvect[k][j]
                        eigenvect[k][j] = s * eigenvect[k][i] + c * eigenvect[k][j]
                        eigenvect[k][i] = vtemp
                    dtemp = c * c * eigenval[i] + s * s * eigenval[j] - 2.0 * c * s * b
                    eigenval[j] = s * s * eigenval[i] + c * c * eigenval[j] + 2.0 * c * s * b
                    eigenval[i] = dtemp
    maxsweeps = m
    for j in range(3):
        k = j
        dtemp = eigenval[k]
        for i in range(j + 1, 4):
            if eigenval[i] < dtemp:
                k = i
                dtemp = eigenval[k]
        if k > j:
            eigenval[k] = eigenval[j]
            eigenval[j] = dtemp
            for i in range(4):
                dtemp = eigenvect[i][k]
                eigenvect[i][k] = eigenvect[i][j]
                eigenvect[i][j] = dtemp
    return eigenvect, eigenval, maxsweeps


def q2mat(quaternion):
    """
    Q2MAT
    Generate a left rotation matrix from a normalized quaternion

    INPUT
      quaternion (q)      - normalized quaternion

    OUTPUT
      rotmat (u)      - the rotation matrix
    """
    rotmat = [[float(0.0) for x in range(3)] for x in range(3)]
    rotmat[0][0] = quaternion[0] * quaternion[0] + quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[1][0] = 2.0 * (quaternion[1] * quaternion[2] - quaternion[0] * quaternion[3])
    rotmat[2][0] = 2.0 * (quaternion[1] * quaternion[3] + quaternion[0] * quaternion[2])
    rotmat[0][1] = 2.0 * (quaternion[2] * quaternion[1] + quaternion[0] * quaternion[3])
    rotmat[1][1] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] + quaternion[2] * quaternion[2] - \
                   quaternion[3] * quaternion[3]
    rotmat[2][1] = 2.0 * (quaternion[2] * quaternion[3] - quaternion[0] * quaternion[1])
    rotmat[0][2] = 2.0 * (quaternion[3] * quaternion[1] - quaternion[0] * quaternion[2])
    rotmat[1][2] = 2.0 * (quaternion[3] * quaternion[2] + quaternion[0] * quaternion[1])
    rotmat[2][2] = quaternion[0] * quaternion[0] - quaternion[1] * quaternion[1] - quaternion[2] * quaternion[2] + \
                   quaternion[3] * quaternion[3]
    return rotmat


def qtrfit(fit_xyz, ref_xyz, maxsweeps):
    """
     QTRFIT
     Find the quaternion, q,[and left rotation matrix, u] that minimizes

       |qTXq - Y| ^ 2  [|uX - Y| ^ 2]

     This is equivalent to maximizing Re (qTXTqY).

     This is equivalent to finding the largest eigenvalue and corresponding
     eigenvector of the matrix

     [A2   AUx  AUy  AUz ]
     [AUx  Ux2  UxUy UzUx]
     [AUy  UxUy Uy2  UyUz]
     [AUz  UzUx UyUz Uz2 ]

     where

       A2   = Xx Yx + Xy Yy + Xz Yz
       Ux2  = Xx Yx - Xy Yy - Xz Yz
       Uy2  = Xy Yy - Xz Yz - Xx Yx
       Uz2  = Xz Yz - Xx Yx - Xy Yy
       AUx  = Xz Yy - Xy Yz
       AUy  = Xx Yz - Xz Yx
       AUz  = Xy Yx - Xx Yy
       UxUy = Xx Yy + Xy Yx
       UyUz = Xy Yz + Xz Yy
       UzUx = Xz Yx + Xx Yz

     The left rotation matrix, u, is obtained from q by

       u = qT1q

     INPUT
       n      - number of points
       fit_xyz (x)      - fitted molecule coordinates
       ref_xyz (y)      - reference molecule coordinates
       weights (w)      - weights

     OUTPUT
       quaternion (q)     - the best-fit quaternion
       rotmat (u)         - the best-fit left rotation matrix
       maxsweeps (nr)     - max number of jacobi sweeps

    """
    # Create variables/lists/matrixes
    matrix = [[float(0.0) for x in range(4)] for x in range(4)]
    quaternion = [float(0.0) for x in range(4)]
    # generate the upper triangle of the quadratic form matrix
    xxyx = float(0.0)
    xxyy = float(0.0)
    xxyz = float(0.0)
    xyyx = float(0.0)
    xyyy = float(0.0)
    xyyz = float(0.0)
    xzyx = float(0.0)
    xzyy = float(0.0)
    xzyz = float(0.0)
    for i, _ in enumerate(fit_xyz):
        xxyx = xxyx + fit_xyz[i][0] * ref_xyz[i][0]
        xxyy = xxyy + fit_xyz[i][0] * ref_xyz[i][1]
        xxyz = xxyz + fit_xyz[i][0] * ref_xyz[i][2]
        xyyx = xyyx + fit_xyz[i][1] * ref_xyz[i][0]
        xyyy = xyyy + fit_xyz[i][1] * ref_xyz[i][1]
        xyyz = xyyz + fit_xyz[i][1] * ref_xyz[i][2]
        xzyx = xzyx + fit_xyz[i][2] * ref_xyz[i][0]
        xzyy = xzyy + fit_xyz[i][2] * ref_xyz[i][1]
        xzyz = xzyz + fit_xyz[i][2] * ref_xyz[i][2]

    matrix[0][0] = xxyx + xyyy + xzyz
    matrix[0][1] = xzyy - xyyz
    matrix[1][1] = xxyx - xyyy - xzyz
    matrix[0][2] = xxyz - xzyx
    matrix[1][2] = xxyy + xyyx
    matrix[2][2] = xyyy - xzyz - xxyx
    matrix[0][3] = xyyx - xxyy
    matrix[1][3] = xzyx + xxyz
    matrix[2][3] = xyyz + xzyy
    matrix[3][3] = xzyz - xxyx - xyyy

    # diagonalize c
    eigenvect, eigenval, maxsweeps = jacobi(matrix, maxsweeps)

    # extract the desired quaternion
    quaternion[0] = eigenvect[0][3]
    quaternion[1] = eigenvect[1][3]
    quaternion[2] = eigenvect[2][3]
    quaternion[3] = eigenvect[3][3]

    # generate the rotation matrix
    rotmat = q2mat(quaternion)

    return quaternion, transpose(rotmat), maxsweeps


def rmsd2(ref_xyz, fit_xyz):
    rms = 0.0
    count = 0
    for co1, co2 in zip(ref_xyz, fit_xyz):
        count += 1
        s1 = (co1[0] - co2[0]) ** 2
        s2 = (co1[1] - co2[1]) ** 2
        s3 = (co1[2] - co2[2]) ** 2
        rms += sum([s1, s2, s3])
    rms = sqrt(rms / count)
    return rms


def quatfitGetMolecule(reffilelol, fitfilelol, pairs):
    """
    ### --- Run the classic quatfit --- ###
    This runs quatfit in a updated form
    reffilelol - reference coordinate in Atom format
    fitfilelol - fit coordinate in molecule format
    pairs - the pairs in a list of lists (i.e. [[2, 3], [3, 13], ...])
    """

    # Create the reference coords based on the pairs
    ref_xyz, fit_xyz = createRefGeom(reffilelol, fitfilelol, pairs)

    # Center the reference coords around 0,0,0
    refcenter, ref_xyz = center(ref_xyz, 1, [float(0), float(0), float(0)])
    fitcenter, fit_xyz = center(fit_xyz, 1, [float(0), float(0), float(0)])

    # fit the specified atom coords of the fit to reference
    quaternion, rotmat, maxsweeps = qtrfit(fit_xyz, ref_xyz, 30)

    # subtract coordinates of the center of fitted atoms of the fitted molecule
    # from all atom coordinates of the fitted molecule
    fitcenter, fitfilelol = center(fitfilelol, 2, fitcenter)

    # rotate the fitted molecule by the rotation matrix u
    fitfilelol = rotmol(fitfilelol, rotmat)

    # same with set of fitted atoms of the fitted molecule
    fit_xyz = rotmol(fit_xyz, rotmat)

    # translate atoms of the fitted molecule to the center
    # of fitted atoms of the reference molecule
    refcenter, fitfilelol = center(fitfilelol, 3, refcenter)

    # same with set of fitted atoms of the fitted molecule
    refcenter, fit_xyz = center(fit_xyz, 3, refcenter)

    # translate fitted atoms of reference molecule to their orig. location
    refcenter, ref_xyz = center(ref_xyz, 3, refcenter)

    rms = 0.0
    s = 0.0
    for i in range(len(pairs)):
        d = 0.0
        for j in range(3):
            if j == 0:
                s = ref_xyz[i].x - fit_xyz[i].x
            elif j == 1:
                s = ref_xyz[i].y - fit_xyz[i].y
            elif j == 2:
                s = ref_xyz[i].z - fit_xyz[i].z
            d += s * s
        rms += d

    rms = sqrt(rms / len(pairs))

    return fitfilelol, rms
