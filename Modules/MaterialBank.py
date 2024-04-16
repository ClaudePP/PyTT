# ----------------------------------- Material Banck ----------------------------------------- #
#
# This module has all the functions and necesary information to load the material data. 
# The material is a class, whith all the material properties being properties of this class. 
#

import sys
import os
import numpy as np
from scipy import constants
from Modules import NecessaryVariables as nv


class Material:

    # When a Object of the class material is created this function is called. 
    # 
    def __init__(self,MaterialFileName):
	
	#
	# To create an object we need as an input the name of the file containing the information 
	# of the material. Follow the examples in the MaterialInfo folder. 
	# This creates a dictionary with the keys being the name of the material property and 
	# and the contents are this values. 
	#
	
        h = open(MaterialFileName)

        d_MatInfo = {}
        cont = 0
        for l in h:
            if cont == 0:
                cont += 1
                continue

            else:
                if len(l.split()) == 0:
                    print("Error: There are blank spaces in Material Properties file.")
                d_MatInfo.update({l.split()[0] : l.split()[1]})

        self.name = d_MatInfo["Name:"]
        self.mpoint = float(d_MatInfo["MeltingPoint:"])             # [K]
        self.rho = float(d_MatInfo["Density:"])                     # [g/cm3] Density
        self.Z = float(d_MatInfo["Z:"])                             # Atomic Number
        self.Am = float(d_MatInfo["Am:"])                           # Atomic Mass
        #self.wfun = float(d_MatInfo["WorkFunction:"])               # [eV] Work function, obsolete, does not work anymore
	
	# The sublimation parameters are not fully needed, one can just have an empty value. 
	#	
        if (d_MatInfo["Sublimation_C1:"] != "-" and d_MatInfo["Sublimation_C2:"] != "-"):
            self.Sublimation_C1 = float(d_MatInfo["Sublimation_C1:"])   # Parameter 1 of sublimation formula. 
            self.Sublimation_C2 = float(d_MatInfo["Sublimation_C2:"])    # Parameter 2 of sublimation formula. 
        else: 
            self.Sublimation_C1 = "-"
            self.Sublimation_C2 = "-"

        #  ------------------- Emissivity ----------------------- #
	#
	# The emissivity can be either a constant value or a string. If the value registered is 
	# a string, the program will look for a file with that same name in the MaterialInflo/ParametersWithTemperature
	# folder. This function loads the values of the parameter in another dictionary. Emissivity dictionary 
	# where the temperature is the key and the values are the emissivity value. 
	# 
	
        EmsInput = d_MatInfo["Emissivity:"]; 
        self.D_Ems = {"Temperature": [], "Parameter": []}
        try: 
            emsT = float(EmsInput)
            self.D_Ems["Temperature"] += [300.0]; self.D_Ems["Parameter"] += [emsT]

        except ValueError:
            pathtoems = os.getcwd()+"/MaterialInfo/ParametersWithTemperature/"+EmsInput+".txt"
            f_ems = open(pathtoems,"r")
            for j,l in enumerate(f_ems,0):
                if j > 3:
                    self.D_Ems["Temperature"] += [float(l.split()[0])]
                    self.D_Ems["Parameter"] += [float(l.split()[1])]

            f_ems.close()

        self.epsT = self.GetEmissivity(np.array([[300]]))   # Default at 300 K

        # ----------------------------------- Specific Heat ----------------------------------------- #
	# 
	# Same as the emissivity
	# units: [J/(g*K)]
	
        CpInput = d_MatInfo["SpecificHeat:"]
        self.D_Cp = {"Temperature": [], "Parameter": []}
        try:
            cpT = float(CpInput)
            self.D_Cp["Temperature"] += [300.0]; self.D_Cp["Parameter"] += [cpT]

        except ValueError:
            pathtocp = os.getcwd()+"/MaterialInfo/ParametersWithTemperature/"+CpInput+".txt"
            f_cp = open(pathtocp,"r")
            for j,l in enumerate(f_cp,0):
                if j > 3:
                    self.D_Cp["Temperature"] += [float(l.split()[0])]
                    self.D_Cp["Parameter"] += [float(l.split()[1])]
            f_cp.close()

        self.CpT = self.GetCp(np.array([[300]]))            # [J/gK] Specific Heat
        
        
        
        #---------------------------------Work function --------------------------
        
        #print("MaterialBank: initialization of work function: ",d_MatInfo["WorkFunction:"],nv.Debug)
        WfInput = d_MatInfo["WorkFunction:"]
        if nv.Debug=="Thermionic":
            print("Debug   MaterialBank:WfInput ",WfInput)
        self.D_Wf = {"Temperature": [], "Parameter": []}
        try:
            Wf = float(WfInput)
            self.D_Wf["Temperature"] += [300.0]; self.D_Wf["Parameter"] += [Wf]

        except ValueError:
            pathtowf = os.getcwd()+"/MaterialInfo/ParametersWithTemperature/"+WfInput+".txt"
            if nv.Debug=="Thermionic":
                print("Debug   MaterialBank:pathtowf ",pathtowf)
            f_wf = open(pathtowf,"r")
            for j,l in enumerate(f_wf,0):
                if j > 3:
                    self.D_Wf["Temperature"] += [float(l.split()[0])]
                    self.D_Wf["Parameter"] += [float(l.split()[1])]
            f_wf.close()
            if nv.Debug=="Thermionic":
                print("Debug   MaterialBank:D_Wf ")
                print(self.D_Wf["Temperature"])
                print(self.D_Wf["Parameter"])
        self.wfun = self.GetWf(np.array([[300]]))            # [eV] work function



        # --------------- Conductivity --------------- #

        Kinput = d_MatInfo["Conductivity:"]
        self.D_k = {"Temperature":[], "Parameter": []}
        try:
            kT = float(Kinput)
            self.D_k["Temperature"] += [300.0]; self.D_k["Parameter"] += [kT]
        
        except ValueError:
            pathtok = os.getcwd()+"/MaterialInfo/ParametersWithTemperature/"+Kinput+".txt"
            f_k = open(pathtok,"r")
            for j,l in enumerate(f_k,0):
                if j > 3:
                    self.D_k["Temperature"] += [float(l.split()[0])]
                    self.D_k["Parameter"] += [float(l.split()[1])]
            f_k.close()

        self.con = self.Getk(np.array([[300]]))             # [W/mK] Conductivity

        # ---------------- Sublimation Enthalpy --------------------- #
        Hinput = d_MatInfo["Sublimation_Heat:"]
        
        if (Hinput == "-"):
            Hinput = 0.0 

        self.D_H = {"Temperature":[], "Parameter": []}
        
        try:
            hT = float(Hinput)
            self.D_H["Temperature"] += [300.0]; self.D_H["Parameter"] += [hT]
        
        except ValueError:
            pathtoh = os.getcwd()+"/MaterialInfo/ParametersWithTemperature/"+Hinput+".txt"
            f_h = open(pathtoh,"r")
            for j,l in enumerate(f_h,0):
                if j > 3: 
                    self.D_H["Temperature"] += [float(l.split()[0])]
                    self.D_H["Parameter"] += [float(l.split()[1])]
            f_h.close()

    
        
        self.HT = self.GetH(np.array([[300]]))              # [J/K] Sublimation enthalpy
        

        self.expcoeff = float(d_MatInfo["ExpansionCoeff:"]) # [1/K] Linear Expansion Coefficient

        h.close()


	#  ----------------------------------------------------------------------------------- #
	
	
	# This function is called during the main loop of the simulations to update the values of 
	# the parameter materials. This is the general function that does the search in a dictionary 
	# for a given temperature. Afterwards there are the various functions feeding this ones the
	# corresponding dictionary per parameter
	
    def GetParameterValue(self,D_Par,T):
        # All the parameter search follows the same logic.                         #
        # ------------------------------------------------------------------------ #
        # If there is only one value, the parameter will be considered constant.   #
        # Otherwise, we will do a linear extrapolation of the parameter values.    #
        #       par(T<300) = par(T=300)                                            #
        #                                                                          #
        #       par(T) = a*T + b    Where a and b are calculated with the values   #
        #       of the parameter  and temperature that are the upper and lower     #
        #       values of T.                                                       #
        #                                                                          #
        #       par(T) with T > par["Temperature"][-1]. We calculate T as an       #
        #       extrapolation of the last two available values.                    #
        # ------------------------------------------------------------------------ #

        #print("MaterialBank/GetParameterValue D_Par:",D_Par,T)  # useless, too much output
        if len(D_Par["Temperature"]) == 1:
            Parameter = D_Par["Parameter"][0]
        else:
            for k in range(0,len(D_Par["Temperature"])):
          
                if (T <= D_Par["Temperature"][k]) and (k == 0):
                    Parameter = D_Par["Parameter"][0]

                    break

                elif (T <= D_Par["Temperature"][k]) and (T < D_Par["Temperature"][-1]):

                    a = (D_Par["Parameter"][k]-D_Par["Parameter"][k-1])/(D_Par["Temperature"][k]-D_Par["Temperature"][k-1])
                    b = D_Par["Parameter"][k]-a*D_Par["Temperature"][k]
                    Parameter = a*T + b
                    
                    break

                elif (T > D_Par["Temperature"][-1] ):

                    a = (D_Par["Parameter"][-1]-D_Par["Parameter"][-2])/(D_Par["Temperature"][-1]-D_Par["Temperature"][-2])
                    b = D_Par["Parameter"][-1]-a*D_Par["Temperature"][-1]
                    Parameter = a*T + b

                    break

        #print(Parameter)
        return Parameter

    def GetEmissivity(self,Temp):
        # Emissivity is a numpy array, with the same dimentions as Temp
        Emissivity = Temp**0
        print("Emissivity: ",Emissivity)
        for i in range(0,len(Emissivity)):
            for j in range(0,len(Emissivity[i])):
                Emissivity[i][j] = self.GetParameterValue(self.D_Ems,Temp[i][j])

        return Emissivity

    def GetCp(self,Temp):
        Cp = Temp**0
        print("MaterialBank Cp: ",Cp)
        for i in range(0,len(Cp)):
            for j in range(0,len(Cp[i])):
                Cp[i][j] = self.GetParameterValue(self.D_Cp,Temp[i][j])
        return Cp


    def GetWf(self,Temp):
        Wf = Temp**0  # create Wf array with the same dimensions as temperature array [[1]]
        print("MaterialBank Wf: ",Wf)
        for i in range(0,len(Wf)):
            for j in range(0,len(Wf[i])):
                Wf[i][j] = self.GetParameterValue(self.D_Wf,Temp[i][j])
        return Wf


    def Getk(self, Temp):
        KK = Temp**0
        for i in range(0,len(KK)):
            for j in range(0,len(KK[i])):
                if Temp[i][j] > self.D_k["Temperature"][-1]:
                    KK[i][j] = self.D_k["Parameter"][-1]
                else:
                    KK[i][j] = self.GetParameterValue(self.D_k,Temp[i][j])
        return KK

    def GetH(self, Temp):
        HH = Temp**0
        for i in range(0,len(HH)):
            for j in range(0,len(HH[i])):
                HH[i][j] = self.GetParameterValue(self.D_H,Temp[i][j])
        return HH


