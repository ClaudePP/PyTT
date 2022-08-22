# -------- Execute simulation ---------------- #
# Note: This simulation will be performed with #
# the input file just created.                 #
# ---------------------------------------------#

import numpy as np
import sys
# This files contain all the files necesary for the temperature evolution simulation.
from Modules import LoadingFileFunctions
from Modules import TargetGeometry
from Modules import CoreSimulationFunctions
from Modules import NecessaryVariables as nv
from Modules import TempPhysicalModels 


# ------------------------------------------------------- #

LoadingFileFunctions.LoadInputFile(nv.RealInputFilename)

# ----------------- Generate Geometry ------------------ #

if nv.DetType != "WIRESCAN":
    TargetGeometry.CreateDetector(nv.DetType)
    nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()
else:

    if nv.BeamType != "Gaussian": 
        print("Wire Scanner Simulations with Non Gaussian beams are still not available .... ")
        sys.exit()
        

# ------------------ Execute Simulations --------------- #

if nv.Flag_Temperature == 1:
    if nv.DetType == "SEM":
        vec = CoreSimulationFunctions.TempEvolSEM()
        nv.V_Time = vec[0]
        nv.V_MaximumTemperature = vec[1]
        nv.M_MaxTemp = vec[2]
        nv.V_Current2 = vec[3]
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
        print("   Options: SEM, FOIL, WIRESCAN")
 
        sys.exit()

# ---------------------------------------------------------#
