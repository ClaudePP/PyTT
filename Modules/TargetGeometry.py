# --------------------- Create detector geometry --------------------------------- #
#
#   This means, we will create some vectors or arrays containing the information of the 
#   location in the x,y plane of our detecotor. The space has been discretized so 
#   this vectors contain the information of all the central locations of the different 
#   divisions of the detector. 
#   This functions also calculate other properties such as the volume of the detectors 
#   and their surface. 
#   
#   In the case of the SEM grids/Wire scanners, the user can specify the plane. 
#   the creation of the geometry is done in the same way, however it is then rotated in case 
#   the vertical plane was requested. 


import numpy as np
import math
from Modules import NecessaryVariables as nv
import sys

def CreateDetector(detectortype):
    if detectortype == 'SEM':
        if nv.SEM_Plane == "Horizontal":
            vec = SEMgridHorizontal(nv.SEM_nWires, nv.SEM_wSep, nv.SEM_wLength, nv.SEM_wWidth, nv.SEM_wRes)
            nv.xvec = np.asanyarray(vec[0])
            nv.yvec = np.asanyarray(vec[1])
            nv.eSup = np.asanyarray(vec[2])
            nv.eVol = np.asanyarray(vec[3])
            nv.IntSurf = np.asanyarray(vec[4])
        elif nv.SEM_Plane == "Vertical":
            vec = SEMgridVertical(nv.SEM_nWires, nv.SEM_wSep, nv.SEM_wLength, nv.SEM_wWidth, nv.SEM_wRes)
            nv.xvec = np.asanyarray(vec[0])
            nv.yvec = np.asanyarray(vec[1])
            nv.eSup = np.asanyarray(vec[2])
            nv.eVol = np.asanyarray(vec[3])
            nv.IntSurf = np.asanyarray(vec[4])

    elif detectortype == 'FOIL':
        vec = FOILDefinition(nv.FOIL_xwidth, nv.FOIL_nx, nv.FOIL_ywidth, nv.FOIL_ny, nv.FOIL_zwidth)
        nv.xvec = np.asanyarray(vec[0])
        nv.yvec = np.asanyarray(vec[1])
        nv.eSup = np.asanyarray(vec[2])
        nv.eVol = np.asanyarray(vec[3])
        nv.IntSurf = np.asanyarray(vec[4])

    elif detectortype == 'WIRESCAN':
        if nv.WIRESCAN_Plane == "Horizontal":
            vec = WIRESCAN_H_Definition(nv.WIRESCAN_wLength, nv.WIRESCAN_wWidth,nv.WIRESCAN_wRes, nv.WIRESCAN_wCposX, nv.WIRESCAN_wCposY)
            nv.xvec = np.asanyarray(vec[0])
            nv.yvec = np.asanyarray(vec[1])
            nv.eSup = np.asanyarray(vec[2])
            nv.eVol = np.asanyarray(vec[3])
            nv.IntSurf = np.asanyarray(vec[4])
        elif nv.WIRESCAN_Plane == "Vertical":
            vec = WIRESCAN_V_Definition(nv.WIRESCAN_wLength, nv.WIRESCAN_wWidth, nv.WIRESCAN_wRes, nv.WIRESCAN_wCposX, nv.WIRESCAN_wCposY)
            nv.xvec = np.asanyarray(vec[0])
            nv.yvec = np.asanyarray(vec[1])
            nv.eSup = np.asanyarray(vec[2])
            nv.eVol = np.asanyarray(vec[3])
            nv.IntSurf = np.asanyarray(vec[4])
        elif nv.WIRESCAN_Plane == "Diagonal":
            vec = WIRESCAN_D_Definition(nv.WIRESCAN_wLength, nv.WIRESCAN_wWidth, nv.WIRESCAN_wRes, nv.WIRESCAN_wCposX, nv.WIRESCAN_wCposY)
            nv.xvec = np.asanyarray(vec[0])
            nv.yvec = np.asanyarray(vec[1])
            nv.eSup = np.asanyarray(vec[2])
            nv.eVol = np.asanyarray(vec[3])
            nv.IntSurf = np.asanyarray(vec[4])
            
        
# ----------------------------------------------------------------------------------------------------------- #

def SEMgridHorizontal(Number_of_wires,wire_separation,wirelength,width,resolution):
    '''
    Sem grid used to make a horizontal profile of the beam. So wires are orientated vertically.
    :param Number_of_wires: Number of wires
    :param wire_separation: Wire separation [m]
    :param wirelength: Wire length [m]
    :param resolution: Wire resolution [m]

    :return: xvec, grid of the x axis.
             yvec, grid of the y axis.
             Sup, surface of each part of the surface.
             Vol, volume of each small part of the surface.
             Isuf, Surface reached by the beam.

    '''

    yvec = np.array(np.arange(-wirelength / 2. + resolution / 2., wirelength / 2. - resolution / 2, resolution))
    xvec = []
    for k in range(0, Number_of_wires + 1):

        if Number_of_wires % 2 == 0:
            # Even number of wires
            if k == Number_of_wires / 2:
                continue
            else:
                xpos = (k - (Number_of_wires / 2.)) * (width + wire_separation)
        else:
            if k == Number_of_wires:
                continue
            else:
                xpos = (k - ((Number_of_wires - 1) / 2.)) * (width + wire_separation)
        xvec += [xpos]

    Sup = math.pi * resolution * width  # [m2] Surface
    Vol = math.pi * resolution * (width / 2) ** 2  # [m3] Volumen
    Isuf = width * resolution
    return np.asanyarray(xvec), yvec, Sup, Vol, Isuf

def SEMgridVertical(Number_of_wires,wire_separation,wirelength,width,resolution):
    '''
    Sem grid used to make a vertical profile of the beam. So wires are orientated vertically.
    :param Number_of_wires: Number of wires
    :param wire_separation: Wire separation [m]
    :param wirelength: Wire length [m]
    :param resolution: Wire resolution [m]

    :return: xvec, grid of the x axis.
             yvec, grid of the y axis.
             Sup, surface of each part of the surface.
             Vol, volume of each small part of the surface.
             Isuf, Surface reached by the beam.

    '''
    xvec = np.array(np.arange(-wirelength / 2. + resolution / 2., wirelength / 2. - resolution / 2, resolution))
    yvec = []

    for k in range(0, Number_of_wires + 1):

        if Number_of_wires % 2 == 0:
            # Even number of wires
            if k == Number_of_wires / 2:
                continue
            else:
                ypos = (k - (Number_of_wires / 2.)) * (width + wire_separation)
        else:
            if k == Number_of_wires:
                continue
            else:
                ypos = (k - ((Number_of_wires - 1) / 2.)) * (width + wire_separation)
        yvec += [ypos]

    Sup = math.pi * resolution * width  # [m2] Surface
    Vol = math.pi * resolution * (width / 2) ** 2  # [m3] Volumen
    Isuf = width * resolution
    return xvec, np.asanyarray(yvec), Sup, Vol, Isuf

def FOILDefinition(width_xaxis, Ncells_xaxis, width_yaxis, Ncells_yaxis, width_zaxis):
    '''
    The foil with be studied as grid. The finner the grid the more precise the calculation. 
    But also more time consuming. 
    width_xaxis: Length in x direction [m].
    Ncells_xaxis: Number of divisions x axis. 
    width_yaxis: Length in y direction [m].
    Ncells_yaxis: Number of divisions y axis.
    width_zaxis: Lenght in z direction [m].
    Note: Z axis is not subdivided. Only valid for "thin" targets. 
    '''

    xvec = np.array(np.arange(-width_xaxis/2., width_xaxis/2., width_xaxis/Ncells_xaxis))
    yvec = np.array(np.arange(-width_yaxis/2., width_yaxis/2., width_yaxis/Ncells_yaxis))

    # This is the surface of a small piece of the grid. 2* as it has two surfaces front and back.
    Sup = 2 * width_xaxis/Ncells_xaxis * width_yaxis/Ncells_yaxis  

    # Volum of a smal piece of the grid. 
    Vol = width_yaxis/Ncells_xaxis*width_yaxis/Ncells_yaxis*width_zaxis

    # Surface exposed to beam. Only front.
    Isurf = width_xaxis/Ncells_xaxis*width_yaxis/Ncells_yaxis

    return xvec, yvec, Sup, Vol, Isurf

def WIRESCAN_H_Definition(wirelength, width, resolution, X0, Y0):
    '''
    This aims to measure the horizontal profile. One wire. X0 fixed Y depends on wlength
    The position in the x axis will be moved as the wire moves. 
    '''
    yvec = np.array(np.arange(-wirelength / 2. + resolution / 2., wirelength / 2. - resolution / 2, resolution))
    xvec = np.array([X0])
    Sup = math.pi * resolution * width  # [m2] Surface
    Vol = math.pi * resolution * (width / 2) ** 2  # [m3] Volumen
    Isuf = width * resolution

    return xvec, yvec, Sup, Vol, Isuf

def WIRESCAN_V_Definition(wirelength, width, resolution, X0, Y0):
    '''
    This aims to measure the Vertical profile. One wire. Y0 fixed X depends on wlength
    The position in the y axes will be updated as the wire moves.
    '''
    yvec = np.array([Y0])
    xvec = np.array(np.arange(-wirelength / 2. + resolution / 2., wirelength / 2. - resolution / 2, resolution))
    Sup = math.pi * resolution * width  # [m2] Surface
    Vol = math.pi * resolution * (width / 2) ** 2  # [m3] Volumen
    Isuf = width * resolution

    return xvec, yvec, Sup, Vol, Isuf

def WIRESCAN_D_Definition(wirelength, width, resolution, X0, Y0):
    print("Diagonally Moving wire scanner is still not ready Sorry!")
    sys.exit()

    return 0.0,0.0,0.0,0.0,0.0,0.0
