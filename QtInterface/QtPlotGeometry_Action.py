########################################################
#                    PLOT GEOMETRY WINDOW              #
#############################################################
#                                                           #
#   In this file the actions in the plot geometry window    #
#   are implemented.                                        #
#                                                           #
#############################################################

# PyQt Libraries
from PyQt5.QtWidgets import  QMainWindow, QApplication, QWidget, QPushButton, QVBoxLayout
from PyQt5.QtCore import *

from QtInterface import QTPlotGeometry_Description

# Python Libraries
import numpy as np
import matplotlib
from matplotlib.figure import Figure

# Necesary libraries for plotting stuff.
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib import pyplot as plt
from matplotlib.patches import Ellipse

# Necessary Things from program itself. 
from Modules import NecessaryVariables as nv


class PlotGeometryWindow(QWidget):
    def __init__(self, InfoListGeo):
        
        super(PlotGeometryWindow, self).__init__()
        self.InfoList = InfoListGeo.copy()
        self.ui = QTPlotGeometry_Description.Ui_Form()
        self.ui.setupUi(self)

        self.figure = Figure()
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout = QVBoxLayout()

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.canvas.ax = self.canvas.figure.add_subplot(111)
        self.ui.widgetGeometryPlot.setLayout(layout)

        self.canvas.ax.set_title("Detector Geometry and Beam", fontsize=10)
        self.canvas.ax.set_xlabel("X [mm]", fontsize = 10)
        self.canvas.ax.set_ylabel("Y [mm]", fontsize = 10)
        X0 = float(self.InfoList[0])*1e3; Y0 = float(self.InfoList[1])*1e3
        Sx = float(self.InfoList[2])*1e3; Sy = float(self.InfoList[3])*1e3
        print(X0, Y0, Sx, Sy)
        self.canvas.ax.add_artist(Ellipse((X0, Y0),Sx,Sy))
        
        self.canvas.draw()


    def WriteGeoText(self):

        self.ui.LabeldettypeGeo.setText(nv.DetType)
        self.ui.LabeldetmatGeo.setText(nv.Material.name)
        self.ui.LabelmtempGeo.setText(nv.Material.mpoint)
        if nv.DetType == "SEM":
            self.ui.labelGeoNumberWires.setText(nv.SEM_nWires)
            self.ui.labelGeoMplane.setText(nv.SEM_Plane)
        elif nv.DetType == "WIRESCAN":
            self.ui.labelGeoNumberWires.setText("1")
            self.ui.labelGeoMplane.setText(nv.WIRESCAN_Plane)
        elif nv.DetType == "FOIL":
            self.ui.labelGeoNumberWires.setText("Foil")
            self.ui.labelGeoMplane.setText("Foil")
        else: 
            self.ui.labelGeoNumberWires.setText(" - ")
            self.ui.labelGeoMplane.setText(" - ")

        
        
    

