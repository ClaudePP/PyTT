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