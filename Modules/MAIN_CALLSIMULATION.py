# ------------- This File launches the simulation --------------------------------- # 
#
# When this module is imported, it starts the simulation.
# If there was an input file given at the beggining, this module was imported in the MAIN.py file. 
# If there was no input file, this module is imported when pressing the simulate button in the user friendly interface. 


# First as always one must load all the necessary modules. 

import numpy as np
import sys

from Modules import LoadingFileFunctions
from Modules import TargetGeometry
from Modules import CoreSimulationFunctions
from Modules import NecessaryVariables as nv
from Modules import TempPhysicalModels 


# ------------- Load Input Files --------------------------------------------------- #
# 
# First of all we load all the necessary parameters given by the input file. 
#
print("MAIN_CALLSIMULATION debug: nv.RealInputFilename is:",nv.RealInputFilename)
LoadingFileFunctions.LoadInputFile(nv.RealInputFilename)

#
# Here We generate the detector geometry. This means:
#    We generate the vectors with the spatial coordinates of the detector. 
#    We calculate some values like the detector volume and surface. 
#    If the detector is a wire scanner this is slightly different, as the coordinates of the detector
#    Will change with the movement. The detector coordinates are created in the simulation loop. 
# The ParticleProportionMatrix calculates for the just calculted detector coordinates, the proportion of beam particles
# hitting that postion. 



if nv.DetType != "WIRESCAN":
    TargetGeometry.CreateDetector(nv.DetType)
    nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()
else:

    if nv.BeamType != "Gaussian": 
        print("Wire Scanner Simulations with Non Gaussian beams are still not available .... ")
        sys.exit()
        
# ------------- Execute Simulations --------------------------------------------------- #
# 
# Here we are ready to execute the simulation. The executed simulation depends on 
# the detector type. In this function we also store the simulation outputs and we 
# write the simulation results in an outputfile. 
# 

if nv.Flag_Temperature == 1:
    print(nv.DetType)
    if nv.DetType == "SEM":
        vec = CoreSimulationFunctions.TempEvolSEM()
        nv.V_Time = vec[0]
        nv.V_MaximumTemperature = vec[1]
        nv.M_MaxTemp = vec[2]
        nv.V_Current2 = vec[3]
        LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)

    elif nv.DetType == "SPLITTER":
        vec = CoreSimulationFunctions.TempEvolSPLITTER()
        nv.V_Time = vec[0]
        nv.V_MaximumTemperature = vec[1]
        #nv.M_MaxTemp = vec[2]
        #nv.V_Current2 = vec[3]   # we should not care about current?
        LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)
            
    elif nv.DetType == "FOIL":
        vec = CoreSimulationFunctions.TimeEvolFOIL()
        nv.V_Time = vec[0]
        nv.V_MaximumTemperature = vec[1]
        nv.M_FancyTemperature = vec[2]
        nv.V_MaxCurrent1 = vec[3]
        nv.V_Current2 = vec[4]
        LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)
            
          
    elif nv.DetType == "WIRESCAN":
        if nv.WIRESCAN_Type == 1:
            vec = CoreSimulationFunctions.TimeEvolWIRESCAN1()
            nv.V_Time = vec[0]
            nv.V_MaximumTemperature = vec[1]
            nv.M_FancyTemperature= vec[2]
            nv.V_Current2 = vec[3]
            nv.V_Pos = vec[4]
            LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)
                
        elif nv.WIRESCAN_Type == 2:
            vec = CoreSimulationFunctions.TimeEvolWIRESCAN2()
            nv.V_Time = vec[0]
            nv.V_MaximumTemperature = vec[1]
            nv.M_FancyTemperature= vec[2]
            nv.V_Current2 = vec[3]
            nv.V_Pos = vec[4]
            LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)
               
        else: 
            print("Select a type of wire simulation: ")
            print("  1: Continuous Beam.")
            print("  2: Pulsated Beam.")
            sys.exit()
    else:
        print("Select a Correct detector Option: ")
        print("   Options: SEM, FOIL, WIRESCAN, SPLITTER")
 
        sys.exit()

# ---------------------------------------------------------#
