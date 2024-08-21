# ----------------------------------- Bethe Bloch formula ----------------------------------------- #
#
# This module has uses infor in ParticleBank and MaterialBank to compute dEdx
# from Bethe-Bloch formula; For high energies and thin targets it gives large
# error due to lack of estimation of delta-electron impact


import numpy as np
from scipy import constants
from Modules import ParticleBank as pb
from Modules import MaterialBank as mb
from Modules import RelativisticFormulas as re


# physics constant:
K = 0.307075 # constant K in [MeV cm mol^-1]
#m_e = 0.511 # [MeV] - from scipy.constants
m_e = constants.physical_constants["electron mass energy equivalent in MeV"][0]
#print(m_e)
I = 0.000173 # Ionisation energy in MeV (of Silicon? verify!) - to Target dataclass? - add to material properties
#mu = 931.494102 # MeV/c^2 - dalton, unified atomic mass unit
mu = constants.physical_constants["atomic mass constant energy equivalent in MeV"][0]

""" maximum energy transferred in one collision in MeV """
""" the same function in DeltaElectrons module? Move to RelativisticFormulas """
def Tmax(partgamma,ionmass):
    '''
    Parameters
    ----------
    partgamma : float
        Relativistic gamma of the projectile.
    ionmass : float
        Ion mass energy equivalent [MeV].

    Returns
    -------
    float
        Maximum energy transfer
    '''    
    #partgamma = lambda partbeta: 1.0/sqrt(1-pow(partbeta,2))
    partbeta = lambda partgamma: np.sqrt(1.0-pow(1.0/partgamma,2))
    return 2*m_e*pow(partbeta(partgamma)*partgamma,2) / (1+2*partgamma*m_e/ionmass+pow(m_e/ionmass,2))

# from: https://pl.wikipedia.org/wiki/Wzór_Bethego-Blocha
def Tmaxv2(partgamma,ionmass):
    '''
    Parameters
    ----------
    partgamma : float
        Relativistic gamma of the projectile.
    ionmass : float
        Ion mass [kg].

    Returns
    -------
    float
        Maximum energy transfer
    '''
    m_ele=constants.electron_mass  # [kg]
    #partbeta = lambda partgamma: np.sqrt(1.0-pow(1.0/partgamma,2))
    partbeta=re.beta_from_gamma(partgamma)
    nominator=2*m_ele*pow(constants.c,2)*pow(partbeta,2)*partgamma**2
    denominat=1+(2*partgamma*m_ele/ionmass)+pow(m_ele/ionmass,2)
    return nominator/denominat



def BetheBloch(particle, bene, material):
    '''
    Computes dEdx from Bethe-Bloch from projectile beta and target Z and A
    '''
    #target.rho = float(d_MatInfo["Density:"])                     # [g/cm3] Density
    #target.Z = float(d_MatInfo["Z:"])                             
    #Am = float(d_MatInfo["Am:"])                                  

    tZ = material.Z                                             # Atomic Number of target   
    tAm = material.Am                                            # Atomic Mass of the target
    tRho = material.rho                                          # [g/cm3] target material density 
    tI = I                                                       # ionization energy of the target, DO!
    pZ = particle.Nprotons

    # what about particle mass???
    # bene=nv.BEnergy            nv.BEnergy                                   # [MeV] - beam particles kinetic energy
    # print(type(bene))
    partgamma = (bene+(pZ*mu))/(pZ*mu)
    partbeta = lambda partgamma: np.sqrt(1.0-pow(1.0/partgamma,2))
    #partgamma = lambda partbeta: 1.0/sqrt(1-pow(partbeta,2))


    C1=K*pow(pZ,2)*tZ/tAm
    return (C1/pow(partbeta(partgamma),2)) * (0.5*np.log(2*m_e*pow(partbeta(partgamma)*partgamma,2)*Tmax(partgamma,pZ*mu)/pow(tI,2))-pow(partbeta(partgamma),2))



