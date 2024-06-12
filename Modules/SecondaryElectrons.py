#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 29 09:00:04 2024

 ------------------ Secondary Electrons -----------------------------------  
 
 This module contains functions dealing with secondary electron properties.


@author: mariusz.sapinski@psi.ch
"""

from Modules import MaterialBank as mb
from Modules import NecessaryVariables as nv
import numpy as np


def ElectronPenetrationLength():
    '''

    Returns
    -------
    Ls : float
        returns effective electron penetration depth in [cm]

    '''        
    Nv = nv.Na*nv.Material.rho/nv.Material.Am   # [atoms/cm3] atomic density 
    # Ls is an effective penetration distance of secondary electrons 
    # 3.68e-17 is a constant with unit of [cm2], see for instance D.Kramer thesis page 29-30
    Ls = 1.0/(3.68e-17*Nv*nv.Material.Z**(1/3.))  # [cm]
    return Ls


def SecondaryElectronYield():
    '''
    
    Returns
    -------
    SEYp : float
        Secondary Electron Yield

    '''
    dEdxel=nv.enemat*1e+6*nv.Material.rho   # [eV/cm] electronic energy loss
    BEnergy_eV = float(nv.BEnergy)*1e6
    pmassamu = 1.00727647
    Ls=ElectronPenetrationLength()
    # units are:
        # 0.01 - [1/eV]
        # Ls - [cm]
        # dEdxel - [eV/cm]
        # 5.4e-6 - [amu/eV]
    SEYp = 0.01*Ls*dEdxel*(1.0+1.0/(1.0+(5.4e-6*BEnergy_eV/pmassamu)))
    return SEYp


# this is in principle for electrons as primary particles
def SecondaryElectronSpectrum(Eks, sigtrue=1.0828, mutrue=1.6636):
    '''
    
    Parameters
    ----------
    Eks : (array of) floats
        electron energies
        
    sigtrue: float
        distribution parameter for true secondary electrons, material-dependent
        default value corresponds to copper covering LHC beam screen
        
    mutru: float
        distribution parameter for true secondary electrons, material-dependent
        default value corresponds to copper covering LHC beam screen
        
    Returns
    -------
    fracs : (array of) floats
        relative intensities for each energy in Eks

    '''
    const= 2.0/(np.sqrt(2*np.pi)*sigtrue*Eks)
    Eksc=[np.log(e-mutrue) for e in Eks]   # negative log argument, something is wrong
    fracs = [c*np.exp(-1*pow(Ec,2)/(2*sigtrue**2)) for c,Ec in zip(const,Eksc)] 
    #const*np.exp(pow(Eksc,2)/(2*psig**2))
    return fracs