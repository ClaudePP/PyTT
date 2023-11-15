############################################################
#                    PyTT MAIN PROGRAM                     #
############################################################

# --------- Importing the necessary Libraries  ----------- #

import sys
import os
from Modules import BetheBloch
from Modules import NecessaryVariables as nv
from Modules import ParticleBank as pb
from Modules import MaterialBank as mb




# The main function checks how many arguments the user gave when executing the code: 
# if python MAIN.py : Then only 1 argument was given. 
# If Python MAIN.py InputFile.txt: Then two arguments were given. 
# 
# In the case of giving no imput file the program calls the user friendly interface. 
# In case there was an inputfile, the MAIN_CALLSIMULATION is directly called. 


def main():

    # launching from command line
    if len(sys.argv) > 1:
        nv.Flag_QtInterface = 0
        nv.RealInputFilename = sys.argv[1]
        from Modules import MAIN_CALLSIMULATION 
    # temporary
    else:
        target = mb.Material("MaterialInfo/Carbon.txt")
        particle = pb.Particle("ParticleInfo/Proton.txt")
        dEdx=BetheBloch.BetheBloch(particle,450000.0,target) # [MeV*cm2/g]
        print(dEdx,dEdx*2*0.0033) # MeV 
        
        
# The code starts here, when the user executes this MAIN.py file. 

if __name__ == '__main__':
    main()
