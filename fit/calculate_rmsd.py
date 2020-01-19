# coding=utf-8
from math import sqrt

from misc import cart_to_frac, frac_to_cart, matrix_minus_vect, matrix_plus_vect

__doc__ = \
    """
    Calculate Root-mean-square deviation (RMSD) of Two Molecules Using Rotation
    ===========================================================================

    Calculate Root-mean-square deviation (RMSD) between structure A and B, in XYZ
    or PDB format, using transformation and rotation. The order of the atoms *must*
    be the same for both structures.

    For more information, usage, example and citation read more at
    https://github.com/charnley/rmsd

    This file was modified by Daniel Kratzert
    """

__version__ = '1.2.7'

import copy

from fit.quatfit import qtrfit, rotmol


def centroid(X):
    """
    Calculate the centroid from a vectorset X.

    https://en.wikipedia.org/wiki/Centroid
    Centroid is the mean position of all the points in all of the coordinate
    directions.

    C = sum(X)/len(X)

    Parameters
    ----------
    X : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    C : float
        centeroid
    >>> centroid([[-1.45300e-02,  1.66590e+00,  1.09660e-01], [-1.46000e-03,  2.68140e-01,  6.35100e-02],[-2.78130e-01, -2.16050e-01,  1.52795e+00]])
    [-0.09804,     0.57266333,  0.56704]
    """
    s = [0.0, 0.0, 0.0]
    num = 0
    for line in X:
        num += 1
        for n, v in enumerate(line):
            s[n] += v
    return (s[0] / num, s[1] / num, s[2] / num)


def rmsd(V, W):
    """
    Calculate Root-mean-square deviation from two sets of vectors V and W.

    Parameters
    ----------
    V : np.array
        (N,D) matrix, where N is points and D is dimension.
    W : np.array
        (N,D) matrix, where N is points and D is dimension.

    Returns
    -------
    fit : float
        Root-mean-square deviation

    """
    D = len(V[0])
    N = len(V)
    rmsd = 0.0
    for v, w in zip(V, W):
        rmsd += sum([(v[i] - w[i]) ** 2.0 for i in range(D)])
    return sqrt(rmsd / N)


def show_coordinates(atoms, V):
    """
    Print coordinates V with corresponding atoms to stdout in XYZ format.

    Parameters
    ----------
    atoms : np.array
        List of atomic types
    V : np.array
        (N,3) matrix of atomic coordinates
    """
    for n, atom in enumerate(atoms):
        atom = atom[0].upper() + atom[1:]
        print("{0:2s}   1 {1:15.8f} {2:15.8f} {3:15.8f}   11.0  0.04".format(atom, *V[n]))


def fit_fragment(fragment_atoms, source_atoms, target_atoms):
    """
    Takes a list of fragment atoms and fits them to the position of target atoms. source_atoms are a fraction
    of the fragment to be fitted on the target_atoms.

    Parameters
    ----------
    fragment_atoms: list
        complete set of atoms of a fragment
    source_atoms: list
        subsection of fragment atoms
    target_atoms: list
        target position for source_atoms

    Returns
    -------
    rotated_fragment: list
        list of coordinates from the fitted fragment
    rmsd: float
        RMSD (root mean square deviation)
    """
    source_atoms = source_atoms
    target_atoms = target_atoms
    fragment_atoms = fragment_atoms
    P_source = copy.deepcopy(source_atoms)
    Q_target = copy.deepcopy(target_atoms)
    # Create the centroid of P_source and Q_target which is the geometric center of a
    # N-dimensional region and translate P_source and Q_target onto that center:
    # http://en.wikipedia.org/wiki/Centroid
    Pcentroid = centroid(P_source)
    Qcentroid = centroid(Q_target)
    # Move P_source and Q_target to origin:
    P_source = matrix_minus_vect(P_source, Pcentroid)
    Q_target = matrix_minus_vect(Q_target, Qcentroid)
    # get the Kabsch rotation matrix:
    quaternion, U, maxsweeps = qtrfit(P_source, Q_target, 30)
    # translate source_atoms onto center:
    source_atoms = matrix_minus_vect(source_atoms, Pcentroid)
    # rotate fragment_atoms (instead of source_atoms):
    rotated_fragment = rotmol(fragment_atoms, U)
    # move fragment back from zero (be aware that the translation is still wrong!):
    rotated_fragment = matrix_plus_vect(rotated_fragment, Qcentroid)
    rms = rmsd(Q_target, P_source)
    return list(rotated_fragment), rms


def test():
    """
    >>> test()
    Kabsch RMSD:   0.0191
    C0d   1      0.13342089      0.24130563      0.55100894   11.0  0.04
    C1d   1      0.17619888      0.17902262      0.56452914   11.0  0.04
    C2d   1      0.08365746      0.12997794      0.52724057   11.0  0.04
    C3d   1     -0.02665530      0.12466866      0.55597929   11.0  0.04
    C4d   1      0.13521359      0.07155518      0.52290434   11.0  0.04
    C5d   1      0.05529952      0.15222482      0.46565980   11.0  0.04
    C6d   1      0.31578577      0.17095916      0.54271372   11.0  0.04
    C7d   1      0.38815031      0.22093469      0.56283111   11.0  0.04
    C8d   1      0.31213535      0.16969488      0.47643097   11.0  0.04
    C9d   1      0.37186868      0.11694406      0.56569282   11.0  0.04
    C10d   1      0.17345616      0.17029947      0.64030933   11.0  0.04
    C11d   1      0.27028472      0.20219330      0.67212698   11.0  0.04
    C12d   1      0.18465947      0.10785452      0.65628617   11.0  0.04
    C13d   1      0.06399788      0.19258933      0.66104906   11.0  0.04

    """
    cell = [10.5086, 20.9035, 20.5072, 90.0, 94.13, 90.0]

    target_atoms = [[0.156860, 0.210330, 0.529750],  # O1
                    [0.198400, 0.149690, 0.543840],  # C1
                    [0.1968, 0.1393, 0.6200]]  # Q11
    # F8: [0.297220,    0.093900,    0.636590]
    # Q11: [0.1968, 0.1393, 0.6200]  # Q11
    target_atoms = [frac_to_cart(x, cell) for x in target_atoms]

    fragment_atoms = [[-0.01453, 1.66590, 0.10966],  # O1  *0
                      [-0.00146, 0.26814, 0.06351],  # C1  *1
                      [-1.13341, -0.23247, -0.90730],  # C2
                      [-2.34661, -0.11273, -0.34544],  # F1
                      [-0.96254, -1.50665, -1.29080],  # F2
                      [-1.12263, 0.55028, -2.01763],  # F3
                      [1.40566, -0.23179, -0.43131],  # C3
                      [2.38529, 0.42340, 0.20561],  # F4
                      [1.53256, 0.03843, -1.75538],  # F5
                      [1.57833, -1.55153, -0.25035],  # F6
                      [-0.27813, -0.21605, 1.52795],  # C4  *10
                      [0.80602, -0.03759, 2.30431],  # F7
                      [-0.58910, -1.52859, 1.53460],  # F8
                      [-1.29323, 0.46963, 2.06735]]  # F9

    target_names = ['O1', 'C1', 'Q11']
    # Structure A rotated and translated onto B (p onto q)

    fragment_atom_names = []
    for n in range(len(fragment_atoms)):
        at = "C{}d".format(n)
        fragment_atom_names.append(at)
    fragment_atom_names = fragment_atom_names
    source_atoms = [fragment_atoms[0], fragment_atoms[1], fragment_atoms[10]]
    rotated_fragment, rmsd = fit_fragment(fragment_atoms, source_atoms, target_atoms)
    rotated_fragment = [cart_to_frac(x, cell) for x in rotated_fragment]  # back to fractional coordinates
    print('Kabsch RMSD: {0:8.3}'.format(rmsd))
    show_coordinates(fragment_atom_names, rotated_fragment)


if __name__ == "__main__":
    test()
