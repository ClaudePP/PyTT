# -*- coding: utf-8 -*-

# This module contains functions to compute ther electric and magnetic field of particle bunch.
# A bunch electric field model must provide the electric field in the REST FRAME of the bunch!
#     .. note::
#       The bunch magnetic field is computed via Lorentz transformation from the bunch
#       electric field (still to be done).
# based on Virtual-IPM code
# (virtual_ipm/simulation/beams/bunches/electric_field_models.py)


import numpy
import scipy.constants as constants
from Modules import NecessaryVariables as nv

"""
    Electrical field for a 2D circular (that is radially symmetric) Gaussian bunch
    (``sigma_x == sigma_y`` within a tolerance of ``%.1e``). <- comment from Virtual_IPM
    here we don't use tolerance
    Solution is obtained from solving Poisson’s equation in 2D.
""" 
def SymmetricBunchField(posr,posz,sigmaz,bcharge):
    """

    Parameters
    ----------
    posr : float
        radial distance from the bunch center
    posz : float
        position along the bunch
    sigmaz : float
        bunch longitudinal sigma (in the bunch frame)
    bcharge : sci float
        number of elementary charges in the bunch

    Returns
    -------
    float
        value of radial electric field in the bunch frame

    """
    sigma=nv.sigx
    inv_sqrt_pi3_eps0_sigmaz = (
        1. / (numpy.sqrt(2. * constants.pi)**3 * constants.epsilon_0 * sigmaz)
    )
    r = posr
    z = posz
    #r_not_equals_zero = (r != 0.) # rm
    # Electric field for r == 0 is zero.
    if posr==0:
        electric_field_abs_val = 0.0
    else:    
        electric_field_abs_val = (
            inv_sqrt_pi3_eps0_sigmaz
            * constants.elementary_charge
            / r
            * (1.0 - numpy.exp(- r**2 / (2. * sigma**2)))
            * numpy.exp(- z**2 / (2. * sigmaz**2))
            )
    #print(electric_field_abs_val)
    return electric_field_abs_val*bcharge


# next:
"""
    Electrical field for a 2D elliptical Gaussian bunch (that is ``sigma_x != sigma_y``). The
    longitudinal charge distribution is taken into account by rescaling the field according to the
    longitudinal position.

    References
    ----------
    .. [1] M.Bassetti & G.A.Erskine: "Closed Expression for the Electrical Field of
       a Two-Dimensional Gaussian Charge", CERN-ISR-TH/80-06, Geneva (CERN), 1980
"""

