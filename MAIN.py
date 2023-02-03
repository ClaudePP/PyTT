############################################################
#                    PyTT MAIN PROGRAM                     #
############################################################

# --------- Importing the necessary Libraries  ----------- #

import sys
import os
from Modules import TargetGeometry 
from Modules import LoadingFileFunctions
from Modules import CoreSimulationFunctions
from Modules import NecessaryVariables as nv

# Set up necessary paths. 
def resource_path(relative_path):
    """ Get absolute path to resource, works for dev and for PyInstaller """
    base_path = getattr(sys, '_MEIPASS', os.path.dirname(os.path.abspath(__file__)))
    return os.path.join(base_path, relative_path)

nv.ParticleInfo_dir = resource_path("ParticleInfo")


# The main function checks how many arguments the user gave when executing the code: 
# if python MAIN.py : Then only 1 argument was given. 
# If Python MAIN.py InputFile.txt: Then two arguments were given. 
# 
# In the case of giving no imput file the program calls the user friendly interface. 
# In case there was an inputfile, the MAIN_CALLSIMULATION is directly called. 


def main():

    if len(sys.argv) > 1:
        nv.Flag_QtInterface = 0
        nv.RealInputFilename = sys.argv[1]
        from Modules import MAIN_CALLSIMULATION 
    else:
        nv.Flag_QtInterface = 1
        from PyQt5.QtWidgets import QApplication, QWidget, QDialog, QMainWindow
        from QtInterface import QtMainWindow_Action
        app = QApplication(sys.argv)
        myapp = QtMainWindow_Action.AppWindow()
        myapp.show()
        sys.exit(app.exec_())
        
        
# The code starts here, when the user executes this MAIN.py file. 

if __name__ == '__main__':
    main()
