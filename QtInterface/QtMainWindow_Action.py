###################################################
#                   MAIN WINDOW                   #
###################################################
#                                                 #
#   This file Opens the main window it contains   #
#   all the functions necesaries for using this   #
#   main window. Those are the button actions.    #
#                                                 #
###################################################

# ------------------ IMPORTING MODULES ------------------------ #

# Imports necessary Python Modules
import sys
import numpy as np
# This imports the necessary PyQt Modules
from PyQt5.QtWidgets import  QMainWindow, QApplication, QWidget, QPushButton
from PyQt5.QtCore import *
# This allows to read html
import codecs

# This imports all the necesary action files that control the different windows. 
from QtInterface import QtMainWindow_Description
from QtInterface import QtPlotGeometry_Action
from QtInterface import QtPlotResultsWindow_Action
from Modules import NecessaryVariables as nv

# -------------- MAIN WINDOW DEFINITION ----------------------# 

class AppWindow(QMainWindow):

    def __init__(self):
        super(AppWindow, self).__init__()
        self.ui = QtMainWindow_Description.Ui_MainWindow()
        self.ui.setupUi(self)

        # This function is executed when Upload File Button is pressed.
        self.ui.pushButtonUploadInputFile.clicked.connect(self.WriteVarToInterface)
        
        # This function is executed when ClearAll button is pressed. 
        self.ui.pushButtonClearAll.clicked.connect(self.DeleteAllText)

        # This actions are executed when the corresponding help buttons are pressed.
        self.ui.pushButtonSUperHelp.clicked.connect(self.GeneralHelp)

        # This action is executed when the PlotGeometry Button is pressed. 
        self.ui.pushButtonPlotGeometrz.clicked.connect(self.PlotGeometryFuncion)


        # This action is executed when the simulate button is pressed. 
        self.ui.pushButtonSimulate.clicked.connect(self.SimulateFunction)


    # ----------------- NECESARY PRESS BUTTON FUNCTIONS ------------------------------ #


    # This function writes automatically the simulation parameters in the user interface.
    def WriteVarToInterface(self):
        wti_Filename = self.ui.lineEditParametersFileName.text()
        f = open(wti_Filename)

        count = 0           # Line counter for FileName
        d_Params = {}       # Dictionary of Necesary Parameters.

        for l in f:
            if len(l.split()) >= 2:
                d_Params.update({l.split()[0] : l.split()[1]})
            else: continue

        f.close()

        # Writes the Beam parameters.
        self.ui.lineEditBeamType.setText(d_Params["BeamType:"])
        if d_Params["BeamType:"] == "Gaussian":
            self.ui.lineEditBeamPosx.setText(d_Params["x0:"])
            self.ui.lineEditBeamPosy.setText(d_Params["y0:"]) 
            self.ui.lineEditBeamSigx.setText(d_Params["sigx:"])
            self.ui.lineEditBeamSigy.setText(d_Params["sigy:"])
        self.ui.lineEditBeamEnergy.setText(d_Params["Energy:"])
        self.ui.lineEditBeamIntensity.setText(d_Params["BCurrent:"])
        self.ui.lineEditNpart.setText(d_Params["Npart:"])
        self.ui.lineEditPlength.setText(d_Params["tpulse:"])
        self.ui.lineEditBeamFrequency.setText(d_Params["frec:"])
        #self.ui.lineEditParticleFile.setText("ParticleInfo/"+d_Params["Particle:"]+".txt") # 2023.10.12: this goes to LoadingFileFunctions
                                                                                            # otherwise command line launching does not work 
        self.ui.lineEditParticleFile.setText(d_Params["Particle:"])


        # Writes Detector Definition Parameters. 

        self.ui.lineEditMaterialFile.setText("MaterialInfo/"+d_Params["Material:"]+".txt")
        self.ui.lineEditMaterialFile.setText(d_Params["Material:"])  # 2023.10.12, moved to LoadingFile

        EneDepFilename = "EneDepData/" + d_Params["Particle:"] + "_" + d_Params["Material:"] + ".txt"
        self.ui.lineEditEnergzDepFile.setText(EneDepFilename)

        if d_Params["DetType:"] == "SEM":
            self.ui.checkBoxSEM.setCheckState(2)
            if d_Params["SEM_Plane:"] == "Horizontal":
                self.ui.checkBoxSEM_H.setCheckState(2)
            elif d_Params["SEM_Plane:"] == "Vertical":
                self.ui.checkBoxSEM_V.setCheckState(2)
            self.ui.lineEditSEMNwires.setText(d_Params["SEM_nWires:"])
            self.ui.lineEditSEMWwidth.setText(d_Params["SEM_wWidth:"])
            self.ui.lineEditSEMWlength.setText(d_Params["SEM_wLength:"])
            self.ui.lineEditSEMWsep.setText(d_Params["SEM_wSep:"])
            self.ui.lineEditSEMWres.setText(d_Params["SEM_wRes:"])
            self.ui.lineEditBeamNpulses.setText(d_Params["SEM_Npulses:"])
        elif d_Params["DetType:"] == "FOIL":
            self.ui.checkBoxFOIL.setCheckState(2)
            self.ui.lineEditFOILxWidth.setText(d_Params["FOIL_xwidth:"])
            self.ui.lineEditFOILnx.setText(d_Params["FOIL_nx:"])
            self.ui.lineEditFOILyWidth.setText(d_Params["FOIL_ywidth:"])
            self.ui.lineEditFOILny.setText(d_Params["FOIL_ny:"])
            self.ui.lineEditFOILzWidth.setText(d_Params["FOIL_zwidth:"])
            self.ui.lineEditBeamNpulses_2.setText(d_Params["FOIL_Npulses:"])
        elif d_Params["DetType:"] == "WIRESCAN":
            self.ui.checkBoxWIRESCANNER.setCheckState(2)
            if d_Params["WIRESCAN_Plane:"] == "Horizontal":
                self.ui.checkBoxWIRE_H.setCheckState(2)
            elif d_Params["WIRESCAN_Plane:"] == "Vertical":
                self.ui.checkBoxWIRE_V.setCheckState(2)

            self.ui.lineEditWIREinipos.setText(d_Params["WIRESCAN_IniPos:"])
            self.ui.lineEditWIREendpos.setText(d_Params["WIRESCAN_EndPos:"])
            self.ui.lineEditWIREwidth.setText(d_Params["WIRESCAN_wWidth:"])
            self.ui.lineEditWIRElength.setText(d_Params["WIRESCAN_wLength:"])
            self.ui.lineEditWIREres.setText(d_Params["WIRESCAN_wRes:"])
            if d_Params["WIRESCAN_Type:"] == "1":
                self.ui.checkBoxWIRE1.setCheckState(2)
                self.ui.lineEditWIREspeed.setText(d_Params["WIRESCAN_wSpeed:"])
            elif d_Params["WIRESCAN_Type:"] == "2":
                self.ui.checkBoxWIRE2.setCheckState(2)
                self.ui.lineEditBeamNpulses_3.setText(d_Params["WIRESCAN_Npulses:"])

        # ---------------------------- Writes Other simulation Parameters. -------------------------------- #
        if d_Params["TempSIM:"] == 'Yes':
            self.ui.checkBoxTEMPERATUREsim.setCheckState(2)
        if d_Params["RadiativeCooling:"] == 'Yes':
            self.ui.checkBoxRadiativeCooling.setCheckState(2)
        if d_Params["ThermionicCooling:"] == "Yes":
            self.ui.checkBoxThermioniCooling.setCheckState(2)
        if d_Params["ConductiveCooling:"] == "Yes":
            self.ui.checkBoxConductionCooling.setCheckState(2)
        if d_Params["SublimationCooling:"] == "Yes":
            self.ui.checkBoxSublimationCooling.setCheckState(2)

        self.ui.lineEditBeamT0.setText(d_Params["T0:"])
        self.ui.lineEditdtPulse.setText(d_Params["dtPulse:"])
        self.ui.lineEditdtCooling.setText(d_Params["dtCooling:"])

        if d_Params["EnableParameterVariation:"] == "Yes":
            self.ui.checkBoxParameterVariation.setCheckState(2)

        if d_Params["IntSIM:"] == "Yes":
            self.ui.checkBoxTEMPERATUREsim_2.setCheckState(2)

        self.ui.lineEditmu.setText(d_Params["mu:"])
        self.ui.lineEditeta.setText(d_Params["eta:"])
        self.ui.lineEditBSp.setText(d_Params["BSp:"])
        self.ui.lineEditBSe.setText(d_Params["BSe:"])

    # This function clears the input parameters in the main window. 
    def DeleteAllText(self):
        self.ui.lineEditParametersFileName.setText("")

        # Deletes Beam Parameters. 
        self.ui.lineEditBeamType.setText("")
        self.ui.lineEditBeamPosx.setText("")
        self.ui.lineEditBeamPosy.setText("")
        self.ui.lineEditBeamSigx.setText("")
        self.ui.lineEditBeamSigy.setText("")
        self.ui.lineEditBeamEnergy.setText("")
        self.ui.lineEditBeamIntensity.setText("")
        self.ui.lineEditNpart.setText("")
        self.ui.lineEditBeamFrequency.setText("")
        self.ui.lineEditBeamNpulses.setText("")
        self.ui.lineEditParticleFile.setText("")
        self.ui.lineEditPlength.setText("")

        # Deletes Detector Definition Parameters. 
        self.ui.lineEditMaterialFile.setText("")
        self.ui.lineEditEnergzDepFile.setText("")

        self.ui.checkBoxSEM.setCheckState(0)
        self.ui.checkBoxSEM_H.setCheckState(0)
        self.ui.checkBoxSEM_V.setCheckState(0)
        self.ui.lineEditSEMNwires.setText("")
        self.ui.lineEditSEMWwidth.setText("")
        self.ui.lineEditSEMWlength.setText("")
        self.ui.lineEditSEMWsep.setText("")
        self.ui.lineEditSEMWres.setText("")
        self.ui.lineEditBeamNpulses.setText("")

        self.ui.checkBoxFOIL.setCheckState(0)

        self.ui.lineEditFOILxWidth.setText("")
        self.ui.lineEditFOILnx.setText("")
        self.ui.lineEditFOILyWidth.setText("")
        self.ui.lineEditFOILny.setText("")
        self.ui.lineEditFOILzWidth.setText("")
        self.ui.lineEditBeamNpulses_2.setText("")

        self.ui.checkBoxWIRESCANNER.setCheckState(0)
        self.ui.checkBoxWIRE_H.setCheckState(0)
        self.ui.checkBoxWIRE_V.setCheckState(0)
    
        self.ui.lineEditWIREinipos.setText("")
        self.ui.lineEditWIREendpos.setText("")
        self.ui.lineEditWIREwidth.setText("")
        self.ui.lineEditWIRElength.setText("")
        self.ui.lineEditWIREres.setText("")

        self.ui.checkBoxWIRE1.setCheckState(0)
        self.ui.lineEditWIREspeed.setText("")
        self.ui.checkBoxWIRE2.setCheckState(0)
        self.ui.lineEditBeamNpulses_3.setText("")

        # Deletes Other simulation Parameters.

        self.ui.checkBoxTEMPERATUREsim.setCheckState(0)
        self.ui.checkBoxTEMPERATUREsim_2.setCheckState(0)

        self.ui.checkBoxRadiativeCooling.setCheckState(0)
        self.ui.checkBoxThermioniCooling.setCheckState(0)
        self.ui.checkBoxConductionCooling.setCheckState(0)
        self.ui.checkBoxSublimationCooling.setCheckState(0)
        self.ui.lineEditBeamT0.setText("")
        self.ui.lineEditdtPulse.setText("")
        self.ui.lineEditdtCooling.setText("")
        self.ui.checkBoxParameterVariation.setCheckState(0)

        self.ui.lineEditmu.setText("")
        self.ui.lineEditeta.setText("")
        self.ui.lineEditBSp.setText("")
        self.ui.lineEditBSe.setText("")

        
    def PlotGeometryFuncion(self):
        # ------ Info List ------- #
        # List with necessary information for drawing detector and beam #
        # Currently it only accepts GAUSSIAN beams                      #
        # Infolist: [X0, Y0, SigmaX, SigmaY, (SEM/FOIL/WIRESCAN), .... ]  #
        # [ ... , [[SEM_H, SEM_V],[], [WIRE_H, WIRE_V, WIRE_45]], ... ] #

        InfoList = []
        InfoList += [self.ui.lineEditBeamPosx.text(), self.ui.lineEditBeamPosy.text()]
        InfoList += [self.ui.lineEditBeamSigx.text(), self.ui.lineEditBeamSigy.text()]
        
        if self.ui.checkBoxSEM.checkState():
            InfoList += ["SEM"]
            if self.ui.checkBoxSEM_H.checkState(): InfoList += ["H"]
            elif self.ui.checkBoxSEM_V.checkState(): InfoList += ["V"]
            else: print("SEM grid: Check plane of measurement")
            InfoList += [self.ui.lineEditSEMNwires.text()]
            InfoList += [self.ui.lineEditSEMWlength.text()]
            InfoList += [self.ui.lineEditSEMWres.text()]
            InfoList += [self.ui.lineEditSEMWsep.text()]
            InfoList += [self.ui.lineEditSEMWwidth.text()]
        elif self.ui.checkBoxWIRESCANNER.checkState():
            InfoList += ["WIRESCAN"]
            if self.ui.checkBoxWIRE_H.checkState(): InfoList += ["H"]
            elif self.ui.checkBoxWIRE_V.checkState(): InfoList += ["V"]
            else: print("WIRESCAN: Check Plane of Measurement")
            InfoList += [self.ui.lineEditWIRElength.text()]
            InfoList += [self.ui.lineEditWIREres.text()]
            InfoList += [self.ui.lineEditWIREwidth.text()]

        elif self.ui.checkBoxFOIL.checkState():
            InfoList += ["FOIL"]
            InfoList += [self.ui.lineEditFOILxWidth.text(), self.ui.lineEditFOILyWidth.text()]
            InfoList += [self.ui.lineEditFOILnx.text(), self.ui.lineEditFOILny.text()]

        else: 
            print(" -------------- Select a detector type! -------------- ")
            sys.exit()

        self.plotgeoWindow = QtPlotGeometry_Action.PlotGeometryWindow(InfoList)
        self.plotgeoWindow.show()
    
    # All this functions open the corresponding help window. 
    def GeneralHelp(self):
        print("Pressing Help")
        f = codecs.open("HelpFolder/PyTT2_MainHelp.html", 'r','utf-8')
    
    def SimulateFunction(self):

        # First a new input file is created with information in main window. 

        rf = open("Output/InputFileUsedForSimulation.txt","w")
        rf.write("# ------------- Input File for PyTT ------------------ #\n")
        rf.write("#     Do NOT change the description names              #\n")
        rf.write("#     There must be at least one space between name    #\n")
        rf.write("#	    and value.					                     #\n")
        rf.write("# ---------------------------------------------------- #\n")
        rf.write("\n\n")
        rf.write("# ------------ Define Beam ------------- #\n")
        rf.write("BeamType: "+self.ui.lineEditBeamType.text()+"\n")
        rf.write("Particle: "+self.ui.lineEditParticleFile.text()+" \n")
        rf.write("Energy: "+self.ui.lineEditBeamEnergy.text()+" \n")
        if self.ui.lineEditBeamType.text() == "Gaussian":
            rf.write("sigx: "+self.ui.lineEditBeamSigx.text()+" \n")
            rf.write("sigy: "+self.ui.lineEditBeamSigy.text()+" \n")
            rf.write("x0: "+self.ui.lineEditBeamPosx.text()+" \n")
            rf.write("y0: "+self.ui.lineEditBeamPosy.text()+" \n")
        rf.write("tpulse: "+self.ui.lineEditPlength.text()+" \n")
        rf.write("BCurrent: "+self.ui.lineEditBeamIntensity.text()+" \n")
        rf.write("Npart: "+self.ui.lineEditNpart.text()+" \n")
        rf.write("frec: "+self.ui.lineEditBeamFrequency.text()+" \n")
        rf.write("# ---------------- Define Detector ----------------- #\n")
        rf.write("Material: "+self.ui.lineEditMaterialFile.text()+" \n")
        rf.write("EneDep: "+self.ui.lineEditEnergzDepFile.text()+"\n")


        if self.ui.checkBoxSEM.checkState():
            rf.write("DetType: SEM \n")
            rf.write("# ------------ Parameters for SEM ------------------ # \n")
            if self.ui.checkBoxSEM_H.checkState(): 
                rf.write("SEM_Plane: Horizontal\n")
            elif self.ui.checkBoxSEM_V.checkState():
                rf.write("SEM_Plane: Vertical\n")
            rf.write("SEM_nWires: "+self.ui.lineEditSEMNwires.text()+" \n")
            rf.write("SEM_wWidth: "+self.ui.lineEditSEMWwidth.text()+" \n")
            rf.write("SEM_wLength: "+self.ui.lineEditSEMWlength.text()+" \n")
            rf.write("SEM_wSep: "+self.ui.lineEditSEMWsep.text()+" \n")
            rf.write("SEM_wRes: "+self.ui.lineEditSEMWres.text()+" \n")
            rf.write("SEM_Npulses: "+self.ui.lineEditBeamNpulses.text()+" \n")

        if self.ui.checkBoxFOIL.checkState():
            rf.write("DetType: FOIL \n")
            rf.write("# -------------- Parameters for FOIL --------------- # \n")
            rf.write("FOIL_xwidth: "+self.ui.lineEditFOILxWidth.text()+" \n")
            rf.write("FOIL_nx: "+self.ui.lineEditFOILny.text()+" \n")
            rf.write("FOIL_ywidth: "+self.ui.lineEditFOILyWidth.text()+" \n")
            rf.write("FOIL_ny: "+self.ui.lineEditFOILny.text()+" \n")
            rf.write("FOIL_zwidth: "+self.ui.lineEditFOILzWidth.text()+" \n")
            rf.write("FOIL_Npulses: "+self.ui.lineEditBeamNpulses_2.text()+" \n")

        if self.ui.checkBoxWIRESCANNER.checkState():
            rf.write("DetType: WIRESCAN \n")
            rf.write("# ------------ Parameters for WIRESCAN ------------- # \n")
            if self.ui.checkBoxWIRE_H.checkState():
                rf.write("WIRESCAN_Plane: Horizontal \n")
            elif self.ui.checkBoxWIRE_V.checkState():
                rf.write("WIRESCAN_Plane: Vertical \n")
            rf.write("WIRESCAN_IniPos: "+self.ui.lineEditWIREinipos.text()+" \n")
            rf.write("WIRESCAN_EndPos: "+self.ui.lineEditWIREendpos.text()+" \n")
            rf.write("WIRESCAN_wWidth: "+self.ui.lineEditWIREwidth.text()+" \n")
            rf.write("WIRESCAN_wLength: "+self.ui.lineEditWIRElength.text()+" \n")
            rf.write("WIRESCAN_wRes: "+self.ui.lineEditWIREres.text()+" \n")
            if self.ui.checkBoxWIRE1.checkState():
                rf.write("WIRESCAN_Type: 1 \n")
                rf.write("WIRESCAN_wSpeed: "+self.ui.lineEditWIREspeed.text()+" \n")
            elif self.ui.checkBoxWIRE2.checkState():
                rf.write("WIRESCAN_Type: 2 \n")
                rf.write("WIRESCAN_Npulses: "+self.ui.lineEditBeamNpulses_3.text()+" \n")
        rf.write("# ------------- Temperature Simulation Parameters -------------- # \n")
        if self.ui.checkBoxTEMPERATUREsim.checkState():
            rf.write("TempSIM: Yes \n")
        else:
            rf.write("TempSIM: No \n")
        rf.write("T0: "+self.ui.lineEditBeamT0.text()+" \n")
        rf.write("dtPulse: "+self.ui.lineEditdtPulse.text()+" \n")
        rf.write("dtCooling: "+self.ui.lineEditdtCooling.text()+" \n")
        
        if self.ui.checkBoxParameterVariation.checkState():
            rf.write("EnableParameterVariation: Yes \n")
        else:
            rf.write("EnableParameterVariation: No \n")
        if self.ui.checkBoxRadiativeCooling.checkState():
            rf.write("RadiativeCooling: Yes \n")
        else:
            rf.write("RadiativeCooling: No \n")
        if self.ui.checkBoxThermioniCooling.checkState():
            rf.write("ThermionicCooling: Yes \n")
        else:
            rf.write("ThermionicCooling: No \n")
        if self.ui.checkBoxConductionCooling.checkState():
            rf.write("ConductiveCooling: Yes \n")
        else:
            rf.write("ConductiveCooling: No \n")
        if self.ui.checkBoxSublimationCooling.checkState():
            rf.write("SublimationCooling: Yes \n")
        else:
            rf.write("SublimationCooling: No \n")

        rf.write("# ------------- Intensity Simulation Parameters -------------- # \n")
        if self.ui.checkBoxTEMPERATUREsim_2.checkState():
            rf.write("IntSIM: Yes \n")
        else:
            rf.write("IntSIM: No \n")

        rf.write("mu: "+self.ui.lineEditmu.text()+" \n")
        rf.write("eta: "+self.ui.lineEditeta.text()+"\n")
        rf.write("BSp: "+self.ui.lineEditBSp.text()+"\n")            
        rf.write("BSe: "+self.ui.lineEditBSp.text()+"\n")    
        rf.write("# --------------------------------------------------- #")

        rf.close()

        # ------------- Executes main siumlation -------- #

        nv.RealInputFilename = "Output/InputFileUsedForSimulation.txt"
        from Modules import MAIN_CALLSIMULATION

        # --------- Open Results Window ------------ #
        
        self.ResultsWindow = QtPlotResultsWindow_Action.PlotResultsWindow()
        self.ResultsWindow.show()

        
    