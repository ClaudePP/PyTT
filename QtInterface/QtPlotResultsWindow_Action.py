########################################################
#                    PLOT RESULTS WINDOW              #
#############################################################
#                                                           #
#   In this file the actions in the plot geometry window    #
#   are implemented.                                        #
#                                                           #
#############################################################


from PyQt5.QtWidgets import  QMainWindow, QApplication, QWidget, QPushButton, QVBoxLayout, QLabel
from PyQt5.QtCore import *
from PyQt5 import QtGui, QtCore

from QtInterface import QTPlotResultsWindow_Description

import numpy as np
import matplotlib
from matplotlib.figure import Figure
from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt4agg import NavigationToolbar2QT as NavigationToolbar

# Necesary libraries for plotting stuff.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import pyplot as plt

from Modules import NecessaryVariables as nv
from Modules import LoadingFileFunctions

class PlotResultsWindow(QWidget):
    def __init__(self):
        super(PlotResultsWindow, self).__init__()
        self.ui = QTPlotResultsWindow_Description.Ui_Form()
        self.ui.setupUi(self)

        # -------- Plot in tab1: Max Temperature vs Time -------- #        
       
        self.PlotResultsMaxTempVStime()

        # -------- Plot in tab2: Intensity vs Time -------- #   

        self.PlotIntensityVsTime()

        # ---------- Plot in tab3: Fancy Temperature  ------------ #

        #self.PlotFancyTemperature()
        #self.PlotMaxMaxTemperature()

        # ----------------------- Add Text ------------------------ #

        self.AddOutputText()

        # ----------------------- Save Button  --------------------- #

        self.ui.pushButtonOutputSave.clicked.connect(self.ConnectWithWritePlots)

        #----------------------------------------------------------- #

        self.canvas.draw()


    def ConnectWithWritePlots(self):
        LoadingFileFunctions.WriteOutputPlotsTxt(nv.OutputFolderName)

    def PlotResultsMaxTempVStime(self):

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.canvas.ax = self.canvas.figure.add_subplot(111)
        self.ui.widgetMaxTempVStime.setLayout(layout)
        self.canvas.ax.grid(color='gray',linestyle="-",linewidth=0.2)
        self.canvas.ax.set_title("Maximum temperature", fontsize=14,fontname="Courier New")
        self.canvas.ax.set_xlabel("Time [s]", fontsize = 14,fontname="Courier New")
        self.canvas.ax.set_ylabel("Temperature [K]", fontsize = 14,fontname="Courier New")
        self.canvas.ax.hlines(nv.Material.mpoint, nv.V_Time[0], nv.V_Time[-1], colors='black', linestyles='solid',lw=1.3)

       
        self.canvas.ax.plot(nv.V_Time,nv.V_MaximumTemperature,color='crimson',lw=1.3)

        #self.canvas.ax.set_ylim([np.min(nv.V_MaximumTemperature[-1])-20.,np.max(nv.V_MaximumTemperature[0])+20.])

    def PlotIntensityVsTime(self):
        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()
        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)
        self.canvas.ax = self.canvas.figure.add_subplot(111)
        self.ui.widgetMaxIntVsTime.setLayout(layout)

        self.canvas.ax.grid(color='gray',linestyle="-",linewidth=0.2)
        self.canvas.ax.set_title("Current in Detector", fontsize=14,fontname="Courier New")
        self.canvas.ax.set_xlabel("Time [s]", fontsize = 14,fontname="Courier New")
        self.canvas.ax.set_ylabel("Current [mA]", fontsize = 14,fontname="Courier New")

        if nv.DetType == "SEM":
            for mk in range(0,len(nv.V_Current2)):
               self.canvas.ax.plot(nv.V_Time,nv.V_Current2[mk]*1e3,color='crimson',lw=1.3) 
        else: 
            self.canvas.ax.plot(nv.V_Time,nv.V_Current2*1e3,color='crimson',lw=1.3)

        


    def PlotFancyTemperature(self):

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout2 = QVBoxLayout()

        layout2.addWidget(self.toolbar)
        layout2.addWidget(self.canvas)

        self.canvas.ax = self.canvas.figure.add_subplot(111)
        self.ui.widgetFanczTemp.setLayout(layout2)

        self.canvas.ax.set_title("Temperature in Detector", fontsize = 16)
        self.canvas.ax.set_xlabel("X [mm]", fontsize = 16)
        self.canvas.ax.set_ylabel("Y [mm]", fontsize = 16)

        if len(nv.xvec) == 1:
            x1 = [nv.xvec[0]-50e-6, nv.xvec[0], nv.xvec[0]+50e-6]
            X,  Y = np.meshgrid(x1, nv.yvec)
            Z = [[0]*len(nv.yvec),list(nv.M_FancyTemperature[0]),[0]*len(nv.yvec)]
        

        elif len(nv.yvec) == 1:
            y1 = [nv.yvec[0]-50e-6, nv.yvec[0],nv.yvec[0]+50e-6]
            X, Y = np.meshgrid(nv.xvec, y1)
            Z = [[0]*len(nv.xvec),list(nv.M_FancyTemperature[0]),[0]*len(nv.xvec)]
        else: 
            X,  Y = np.meshgrid(nv.xvec, nv.yvec)
            Z = nv.M_FancyTemperature

        
        
        cs = self.canvas.ax.contourf(X*1e+3,Y*1e+3,np.transpose(Z), 10)
        
        img = self.canvas.ax.imshow(np.transpose(Z), cmap = "inferno", interpolation='none')
        self.figure.colorbar(img)



    def AddOutputText(self):

            self.ui.Labeldettype.setText(nv.DetType)
            self.ui.Labeldetmat.setText(nv.Material.name)
            self.ui.Labelmtemp.setText(str(nv.Material.mpoint)+" [K]")
            self.ui.labelmaxtemp.setText(str(round(np.max(nv.V_MaximumTemperature),3))+" [K]")
            #self.ui.labelmaxel.setText(str(round(np.max(nv.WireExp)*1e+3,5))+" [mm]")  
