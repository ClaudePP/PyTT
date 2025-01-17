#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 12 08:33:00 2024

@author: mariusz.sapinski@psi.ch
"""

from Modules import NecessaryVariables as nv
import numpy as np
from Modules import MaterialBank as mb
from scipy.constants import physical_constants


# maximum kinetic energy of delta electrons
def Edmax(beta,Eshell=0):
    '''

    Parameters
    ----------
    beta : float
        Relativistic beta of projectile particles (dimensionless, 0 ≤ beta < 1).
    Eshell: float    
        Binding energy of the shell from which the electron is extracted, in [eV]. Default is 0.
    Returns
    -------
    float
        Maximum energy of delta electrons in Joules [J = kg*m^2/s^2]
    Notes
    -----
    - Requires constants: `nv.Qe` (elementary charge in Coulomb), 
      `nv.Emass` (electron mass in kg), and `nv.cspeed` (speed of light in m/s).
    """

    '''
    if not (0 <= beta < 1):
        raise ValueError("Beta must be between 0 (inclusive) and 1 (exclusive).")    
    
    electron_mass = physical_constants["electron mass"][0]  # kg
    speed_of_light = physical_constants["speed of light in vacuum"][0]  # m/s
    elementary_charge = physical_constants["elementary charge"][0]  # C (J/eV)    
    
    
    # Convert binding energy from eV to J    
    #Eshell_J=Eshell*nv.Qe  # [eV]->[J]
    Eshell_J=Eshell*elementary_charge  # [eV]->[J]
    
    
    # Calculate maximum kinetic energy
    max_energy = 2*electron_mass*pow(speed_of_light*beta,2)/(1.0-beta**2)+Eshell_J 

    return max_energy




# Kanaya-Okayama Formula
# for comparison of various formulas see:
# https://www.researchgate.net/publication/225329252_Depth_range_of_primary_electrons_electron_beam_broadening_and_spatial_resolution_in_electron-beam_studies    
def PenetrationDepth(Ee,Zt,At,rho):
    '''
    Parameters
    ----------
    Ee : (array of) floats
        electron energy [eV]
    Zt : integer    
        nuclear charge of target material
    At : float
        atomic number
    rho : float
        material density [g/cm3]    
    Returns
    -------
    Rele : (array of) floats
        electron range in the material [m]

    '''
    Ee=Ee/1000  # convert eV to keV
    Rele=27.6e-9*At*pow(Ee,1.67)
    #if Ee<5.0:
    #    Rele=27.6e-9*At*pow(Ee,1.35)        
    Rele=Rele/(rho*pow(Zt,0.89))
    return Rele




# see: https://pdf.sciencedirectassets.com/271559/1-s2.0-S0168583X13X00308/1-s2.0-S0168583X13011634/main.pdf?X-Amz-Security-Token=IQoJb3JpZ2luX2VjENb%2F%2F%2F%2F%2F%2F%2F%2F%2F%2FwEaCXVzLWVhc3QtMSJHMEUCIBJccF6Y7f8jA3o23dHcdynzfv9yPEjGSy5cMnwMG9NQAiEAsVsmVO1sSOC26AHWixSVGrhKrK1vyrw39N9gdeiTv8oqswUIbxAFGgwwNTkwMDM1NDY4NjUiDE5nKRLzbOwoflH7oCqQBdK78oSrVp%2Fb1uqP9zqBiMLxz7mqhU9LIrO7x2rkgMLjZtZZNiJjQGzpwJRErGN5dDzXjvxeLutsRkbFlZ6wQUYF7q%2FXqJQKVWOiRLUs46aZaztWWYgcZJDwnxxO6N5nLufDrOI3Cl1sG8n%2BxTjQtnpQ6JJFUZkGWQQChox1bIoQcFzCuoAF21IyAHJExC0d02GifJcpTnHnaCqJgBr67QgnVjjXRge5xvJrzhKp5fux0wLTi7RKO93vqOCn0C6Nlc5doMYe4JP2vOY4k7MmMrHTlAJxdvKN81Skvi0M32UoNNPNw8qCoBcGdeb0NDQLGIKBgGksS1jCSPE08JGmmgDsmUM%2BPvLQlvdJQZohddP8iaXFmGAxCVeqRUvYCym6c0%2FwW5kH7OdBdnvjN79ELTd%2FfV86uqN%2FjPuzrDvfSioDN1MeuCJMRmG8JBYtNXX1K7WPWLzZztlJW2%2BkzCcLtcfZFwjZLJKBCYABztAAmSFLmLkzfo8%2FDrNo2fHr5f7OZbr9xcr0iuwNU7aY%2FTfIMMNKMfF9Z%2BBRNsx5qqhXbeCUK27x5bt0tLe4vdQfLc9ZysFna7LEKdU1xE6ANYhwnFhWi%2B2XOdnsyl2hdec7IegKI0IUeLREF0mMmmAyUJpsR8wN3miMSeaqe%2BkcGZARSSaLct8qLNj7xUVuGoTLKKMGODGe%2FcHL%2B5dArKw9tWI8uxkBWkSSR2vDRo3r%2BSNSVEy%2BIbuz98HY8fh9oiESw%2FfPzQSx1uhEUzL%2B4%2BfJ%2F%2BtIxb%2FI2fjcVeIwccxDud9fwQGD0sXKofKZoJLDHrKmaTnGXpD2vq10%2FbdQMtsHw16QncDjbYcUlW3qmr2sDW75bYm4nF88YdSNO%2FvD6WQFHOHZMJ7lpLMGOrEBp6SPNqJcr71qiNm7TWkmLk1crN4EK3xAc6QjewtaUQOToVZR8TOTC4cYMTJVMGnktx3SpiZIrq4NT0doA3psYZTkwt7UKkxlhy9bpopnUHYdjpUxoy1nWdlGmyTWGbjtwoDwmR12I6YcmdAPJInnEkrzHO%2BrtiHYlXFitBIBk95iSs7mk4WTeJ7T8r8Vs4m3E80tu1usupFcM%2B5FxE%2B0aLkNgRWsZHYOcUt%2Fnh0fbKlj&X-Amz-Algorithm=AWS4-HMAC-SHA256&X-Amz-Date=20240612T064526Z&X-Amz-SignedHeaders=host&X-Amz-Expires=300&X-Amz-Credential=ASIAQ3PHCVTYVBFDVPOA%2F20240612%2Fus-east-1%2Fs3%2Faws4_request&X-Amz-Signature=297f22dfdbf4b52cdb421c9601bc993ecb6da5c90d869590d44b0d080360e5d3&hash=184f0da793fa27e8367b5fac1b7d8fb84c7fd13f78b44b80d8d89edd3359a817&host=68042c943591013ac2b2430a89b270f6af2c76d8dfd086a07176afe7c76c2c61&pii=S0168583X13011634&tid=spdf-eefc4b65-d5ed-402e-8a5f-349b27a3d4c1&sid=7c8bdc1d318e5241522b74579bb2a7c7dd5cgxrqb&type=client&tsoh=d3d3LnNjaWVuY2VkaXJlY3QuY29t&ua=05085f5e525803510350&rr=8927e0c7bf566aa1&cc=ch
# A. Akkerman et al. / Nuclear Instruments and Methods in Physics Research B 321 (2014) 1–7
# formula used by Geant4
def DiffECrossSection(Edelta,Zt,beta,Zeff,Eshell=0.0):
    '''

    Parameters
    ----------
    Ee : (array of) floats
        electron energies [eV]
    Zt : integer
        nuclear charge of target atoms
    beta : float
        v/c of the projectile
    Zeff : float
        effective charge of the projectile
    Eshell: float
        [eV] electron binding energy


    Returns
    -------
    float
        (array of) differential cross sections at energies defined by Ee
        in units of [m**2/J]
    '''
    
    Edelta = Edelta*nv.Qe        # convert eV to Joules
    Eshell = Eshell*nv.Qe        # convert eV to Joules
    Ee=Edelta+Eshell
    comp1=2*np.pi*nv.re**2*nv.Emass*nv.cspeed**2*Zt*Zeff**2/pow(beta*Ee,2)
    Emax=Edmax(beta,Eshell)
    # last term in the following is for spin 1/2 (proton OK)
    comp2=1.0-(beta**2*Ee/Emax)+(Ee**2/(2*Emax**2))
    #comp2=1.0-(beta**2*Ee/Emax)
    print("comp1, comp2, comp1*comp2 = ", comp1,comp2,comp1*comp2)
    #print("comp2 = ", comp2)
    #print("Ee = ",Ee)
    #print("Emax = ",Emax)

    return comp1*comp2
    
