#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:44:42 2024

@author: sapinski
"""

import numpy

# copied from VirtualIPM code (author is D.Vilsmeier)
def beta_gamma_from_energy(kin_energy, rest_energy):
    """
    Compute relativistic beta and gamma values from particle energy.
    Both energies have to have the same units.
    
    Parameters
    ----------
    energy : float
        Kinetic energy of the particle.
    rest_energy : float
        Rest energy of the particle.
        
    Returns
    -------
    tuple
        Tuple containing (beta, gamma).
    """
    gamma = 1. + kin_energy / rest_energy
    beta = numpy.sqrt(1. - 1. / gamma ** 2)
    return beta, gamma


def beta_from_gamma(gamma):
    """
    Parameters
    ----------
    gamma : float
        Relativistic gamma.

    Returns
    -------
    float
        Relativistic beta.

    """
    return numpy.sqrt(1.0-pow(1.0/gamma,2))