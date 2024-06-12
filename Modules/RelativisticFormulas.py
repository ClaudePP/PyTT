#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 10:44:42 2024

@author: sapinski
"""

import numpy

def compute_beta_and_gamma_from_energy(energy, rest_energy):
    """
    Compute relativistic beta and gamma values from particle energy.
    
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
    gamma = 1. + energy / rest_energy
    beta = numpy.sqrt(1. - 1. / gamma ** 2)
    return beta, gamma
