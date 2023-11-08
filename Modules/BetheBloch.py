# ----------------------------------- Bethe Bloch formula ----------------------------------------- #
#
# This module has uses infor in ParticleBank and MaterialBank to compute dEdx
# from Bethe-Bloch formula; For high energies and thin targets it gives large
# error due to lack of estimation of delta-electron impact


import numpy as np
from scipy import constants

# physics constant:
K = 0.307075 # constant K in [MeV cm mol^-1]
#m_e = 0.511 # [MeV] - from scipy.constants
m_e = constants["electron mass enery equivalent in MeV"][0]
I = 0.000173 # Ionisation energy in MeV (of Silicon? verify!) - to Target dataclass? - add to material properties
#mu = 931.494102 # MeV/c^2 - dalton, unified atomic mass unit
mu = constants["atomic mass constant energy equivanlent in MeV"]

""" maximum energy transferred in one collision in MeV """
def Tmax(partgamma,ionmass):
    #partgamma = lambda partbeta: 1.0/sqrt(1-pow(partbeta,2))
    partbeta = lambda partgamma: sqrt(1.0-pow(1.0/partgamma,2))
    return 2*m_e*pow(partbeta*partgamma,2) / (1+2*partgamma*m_e/ionmass+pow(m_e/ionmass,2))


def BetheBloch(particle, target, bene)
    '''
    Computes dEdx from Bethe-Bloch from projectile beta and target Z and A
    '''
    target.rho = float(d_MatInfo["Density:"])                     # [g/cm3] Density
    target.Z = float(d_MatInfo["Z:"])                             # Atomic Number
    Am = float(d_MatInfo["Am:"])                                  # Atomic Mass
    bene=nv.BEnergy                                               # [MeV] - beam particles kinetic energy
    partgamma = (bene+(Am*mu))/(Am*mu)
    partbeta = lambda partgamma: sqrt(1.0-pow(1.0/partgamma,2))
    #partgamma = lambda partbeta: 1.0/sqrt(1-pow(partbeta,2))


    C1=K*pow(particle.Z,2)*target.Z/target.Am
    return (C1/pow(partbeta,2)) * (0.5*log(2*m_e*pow(partbeta*partgamma,2)*Tmax(partgamma,ion)/pow(target.I,2))-pow(partbeta,2))



