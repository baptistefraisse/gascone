""" Angular momentum studies functions """


# librairies

from env import *
import numpy as np


# functions 


def angmom_capture(e):

    """
    Compute the evolution of gamma-rays multiplicity from 
    its value at the lowest incident neutron energy.

    Args:
        a (int): Mass number of target nucleus.
        e (float or narray): Incident neutron energy.

    Returns:
        float or narray: Angular momentum from neutron capture  [hbar unit].
    """
    j0 = np.sqrt(2.5*e + SN*SN)
    de = 0.5 # 1 MeV binning
    j0_err = de * 2.5 / (2*j0)
    return j0, j0_err


def angmom_emission(s, nubar):

    """
    Angular momentum take away by neutron emissions.

    Args:
        s (float or narray): Average angular momentum taken away by a single neutron [hbar unit].
        nubar (float or narray): Average neutron multiplicity.

    Returns:
        float or narray: Average angular momentum taken away by the neutron cascade [hbar unit].
    """
    return (s-SN) * nubar


def angmom_frag(frag_ratio, l):
    """
    Compute the angular momentum transmitted to the fragments.

    Args:
        frag_ratio (float): Ratio of angular momentum transmitted to the fragments (0 to 1).
        l (float or narray): Fissionning system angular momentum [hbar unit].

    Returns:
        float or narray: Fission fragments angular momentum [hbar unit].
    """
    return frag_ratio * l


def electric_trans(pole, l):
    """
    Gamma-rays multiplicity of a purely electric decay from fission fragments angular momentum.

    Args:
        pole (int): Electric transition multipole order (1, 2, 3 ...) 
        l (float or narray): Nucleus angular momentum [hbar unit].

    Returns:
        float or narray: Gamma-rays multiplicity cascade.
    """

    return (1./pole) * l


def g_mult_electrans(j0, s, frag_ratio, pole):
    """
    Gamma-rays multiplicity of a purely electric decay from incident neutron angular momentum.

    Args:
        j0 (narray): Initial angular momentum [hbar unit].
        s (float or narray): Average angular momentum taken away by a single neutron [hbar unit].
        frag_ratio (float): Orbital ration of angular momentum transmissions (0 to 1).
        pole (int): Electric transition multipole order (1, 2, 3, ...).

    Returns:
        float or narray: Incident neutron angular momentum.
        float or narray: Incident neutron angular momentum 1-sigma error.
        float or narray: Gamma-rays multiplicity cascade.
    """

    # transmission to the fragments

    j_frag = angmom_frag(frag_ratio, j0)

    # neutron emissions subtraction

    jg = j_frag - angmom_emission(s, nubar_jeff[1:])

    # electric transisions of remnant momentum
    
    g_mult = electric_trans(pole, jg)

    return g_mult