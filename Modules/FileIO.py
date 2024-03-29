# ------------- This File reads the input file and writes the output file --------------------------------- # 

from Modules import NecessaryVariables as nv
from Modules import ParticleBank as pb
from Modules import MaterialBank as mb
from Modules import BetheBloch as bb
import datetime
import sys
import numpy as np

# -------------     Load Input File    ------------------------------------------------------------- #
# 
# To read the input file we will read all the lines in the given file and store them in a dictionary. 
# It is important that the input file follows the format Property: Value
# Including the : and the space. Otherwise this function will not work properly. 
# It is also important that the property names in the input file match the ones expected by this function. 
# Look at the InputFile_Example.txt for a reference

def LoadInputFile(FileName):

    f = open(FileName)
    print("LoadingFunctions debug: FileName is:",FileName)

    count = 0         # Line counter for FileName
    d_Params = {}     # Dictionary of Necesary Parameters.

    # First we read all the lines, ignoring empty ones and save them in the dictionary
    
    for l in f: 
        
        if len(l.split()) >= 2:
            d_Params.update({l.split()[0] : l.split()[1]})
        else: continue
        
        count += 1

    f.close()
    
    #
    # Now we match the values of this dictionary with the expected properties names and 
    # we update the nevesary variables file with this inputs. 
    #
    
    #  1)  Beam Parameters

    nv.BeamType = d_Params["BeamType:"]
    print("LoadingFileFunctions debug: d_Params[Particle:] is ",d_Params["Particle:"])
    nv.Particle = pb.Particle("ParticleInfo/"+d_Params["Particle:"]+".txt")  # here
    nv.BEnergy = d_Params["Energy:"]
    if nv.BeamType == "Gaussian":
        nv.sigx = float(d_Params["sigx:"])
        nv.sigy = float(d_Params["sigy:"])
        nv.x0 = float(d_Params["x0:"])
        nv.y0 = float(d_Params["y0:"])
    nv.tpulse = float(d_Params["tpulse:"])     # [s] beam pulse duration

    if (d_Params["Npart:"] == "-") and (d_Params["BCurrent:"] == "-"):
        print("Beam Current not correctly set: No Value"); sys.exit()
    elif (d_Params["Npart:"] != "-") and (d_Params["BCurrent:"] != "-"):
        print("Beam Current not correctly set: One must give either intensity or number of particles. Not both."); 
        sys.exit()
    elif (d_Params["Npart:"] != "-") and (d_Params["BCurrent:"] == "-"):
        nv.Nparticles = float(d_Params["Npart:"])
        nv.Intensity = 0.0
    elif (d_Params["Npart:"] == "-") and (d_Params["BCurrent:"] != "-"):
        nv.Nparticles = 0.0
        nv.Intensity = float(d_Params["BCurrent:"])

    nv.frec = float(d_Params["frec:"])

    # 2)  Material and Detector Parameters 
    # Notice that only the necessary variables will be read. Even if 
    # the input file contains more lines its ok. 
    
    nv.Material = mb.Material("MaterialInfo/"+d_Params["Material:"]+".txt")
    nv.DetType = d_Params["DetType:"]

    if nv.DetType == "SEM":
        nv.SEM_Plane = d_Params["SEM_Plane:"]
        nv.SEM_nWires = int(d_Params["SEM_nWires:"])
        nv.SEM_wWidth = float(d_Params["SEM_wWidth:"])
        nv.SEM_wLength = float(d_Params["SEM_wLength:"])
        nv.SEM_wSep = float(d_Params["SEM_wSep:"])
        nv.SEM_wRes = float(d_Params["SEM_wRes:"])
        nv.Npulses = int(d_Params["SEM_Npulses:"])
    elif nv.DetType == "SPLITTER":
        nv.SPLITTER_Plane = d_Params["SPLITTER_Plane:"]
        nv.SPLITTER_wPos = d_Params["SPLITTER_wPos:"]
        nv.SPLITTER_wWidth = float(d_Params["SPLITTER_wWidth:"])
        nv.SPLITTER_wDepth = float(d_Params["SPLITTER_wDepth:"])
        nv.SPLITTER_wLength = float(d_Params["SPLITTER_wLength:"])
        nv.SPLITTER_wRes = float(d_Params["SPLITTER_wRes:"])
    elif nv.DetType == "FOIL":
        nv.FOIL_xwidth = float(d_Params["FOIL_xwidth:"])
        nv.FOIL_nx = float(d_Params["FOIL_nx:"])
        nv.FOIL_ywidth = float(d_Params["FOIL_ywidth:"])
        nv.FOIL_ny = float(d_Params["FOIL_ny:"])
        nv.FOIL_zwidth = float(d_Params["FOIL_zwidth:"])
        nv.Npulses = int(d_Params["FOIL_Npulses:"])
    elif nv.DetType == "WIRESCAN":
        nv.WIRESCAN_Plane = d_Params["WIRESCAN_Plane:"]
        nv.WIRESCAN_IniPos = float(d_Params["WIRESCAN_IniPos:"])
        nv.WIRESCAN_EndPos = float(d_Params["WIRESCAN_EndPos:"])
        nv.WIRESCAN_wShape = d_Params["WIRESCAN_wShape:"]
        nv.WIRESCAN_wWidth = float(d_Params["WIRESCAN_wWidth:"])
        if nv.WIRESCAN_wShape == "Strip":                                 # 2023.11.02: this is not yet implemented
            nv.WIRESCAN_wDepth = float(d_Params["WIRESCAN_wDepth"])
        nv.WIRESCAN_wLength = float(d_Params["WIRESCAN_wLength:"])
        nv.WIRESCAN_wRes = float(d_Params["WIRESCAN_wRes:"])
        nv.WIRESCAN_Type = int(d_Params["WIRESCAN_Type:"])
        if nv.WIRESCAN_Type == 1:
            nv.WIRESCAN_wSpeed = float(d_Params["WIRESCAN_wSpeed:"])
        elif nv.WIRESCAN_Type == 2: 
            nv.Npulses = int(d_Params["WIRESCAN_Npulses:"])
        else: 
            print("Wire Scanner Type not selected")
            print("Exiting ...... ")
            sys.exit()
        
    else:
        print("The detector type selected is not available")
        print("Exiting .....")
        sys.exit()
    
    #
    # 3) Temperature Simulation Parameters 
    #
    print("LoadingFileFunctions debug: d_Params is:", d_Params)
    #if d_Params["TemperatureSimulation:"] == 'Yes':
    if d_Params["TempSIM:"] == 'Yes':        
        nv.Flag_Temperature = 1
    #elif d_Params["TemperatureSimulation:"] == 'No':
    else:    
        nv.Flag_Temperature = 0

    nv.T0 = float(d_Params["T0:"])
    nv.dtPulse = float(d_Params["dtPulse:"])           # [s] simulation resolution during the beam pulse
    nv.dtCooling = float(d_Params["dtCooling:"])       # [s] simulation resolution after the beam pulse (cooling only)

    if d_Params["EnableParameterVariation:"] == 'Yes':
        nv.EnableParameterVariation = 1
    elif d_Params["EnableParameterVariation:"] == 'No':
        nv.EnableParameterVariation = 0
    
    if d_Params["RadiativeCooling:"] == 'Yes':
        nv.RadiativeCooling = 1
    elif d_Params["RadiativeCooling:"] == 'No': 
        nv.RadiativeCooling = 0 

    if d_Params["ThermionicCooling:"] == 'Yes':
        nv.ThermoionicCooling = 1
    elif d_Params["ThermionicCooling:"] == 'No': 
        nv.ThermoionicCooling = 0
    
    if d_Params["ConductiveCooling:"] == 'Yes':
        nv.ConductiveCooling = 1
    elif d_Params["ConductiveCooling:"] == 'No':
        nv.ConductiveCooling = 0

    if d_Params["SublimationCooling:"] == 'Yes':
        if (nv.Material.Sublimation_C1 != '-') and (nv.Material.Sublimation_C2 != '-'):
            nv.SublimationCooling = 1
        else: 
            nv.SublimationCooling = 0
    elif d_Params["SublimationCooling:"] == 'No':
        nv.SublimationCooling = 0

    #    
    # 4) Intensity Simulation Parameters 
    #
    
    #if d_Params["IntensitySimulation:"] == 'Yes':
    if d_Params["IntSIM:"] == 'Yes':
        nv.Flag_Intensity = 1
    #elif d_Params["IntensitySimulation:"] == 'No':
    else:
        nv.Flag_Intensity = 0

    nv.Mu = float(d_Params["mu:"])
    nv.Eta = float(d_Params["eta:"])
    nv.BEp = float(d_Params["BSp:"]) 
    nv.BEe = float(d_Params["BSe:"]) 

    #
    # 5)  Load appropiate energy deposition 
    #
    # Notice that the energy deposition given in the input file 
    # is a file name, with its corresponding path. 
    # Here we read this file. We also create a dictionary. In this case 
    # the keys of the dictionary are strings with the energy of the incident particles
    # The nv.BEnergy has to match one of these keys or the program will complain. 
    #
    # See EneDepData folder for examples of Energy deposition files. 
    
    print(d_Params)
    nv.EdepMethod = d_Params["EdepMethod:"]
    if nv.EdepMethod == "Interpolated":
        # data with dE/dx:
        EneDep_Filename = "EneDepData/" + d_Params["Particle:"] + "_" + d_Params["Material:"] + ".txt"
    
        #g = open(EneDep_Filename)
        #d_EneDep = {}   # Dictionary with energy loss 
        #cont = 0
        #for l in g:
        #    if cont == 0:
        #        cont+=1
        #        continue
        #    else:
        #        d_EneDep.update({l.split()[0] : l.split()[1]})        
        #nv.enemat = float(d_EneDep[nv.BEnergy])  # do the interpolation!
        #g.close()
        beamef=[]
        dedxf=[]
        with open(EneDep_Filename) as dedxfile:
            for idx,line in enumerate(dedxfile):
                print(idx,line)
                if idx>0:
                    beame,dedx = line.split()
                    beamef.append(float(beame))
                    dedxf.append(float(dedx))
        #print(beamef,dedxf)            
        # check for exceptions             
        nv.enemat = np.interp(nv.BEnergy,beamef,dedxf)
    elif nv.EdepMethod == "Value":
        # first check if field EneDep exists, throw exception if not
        try: 
            nv.enemat = float(d_Params["Edep:"])
        except:
            print("Missing Edep field describing dE/dx [MeV*cm2/g] for this beam.")
            exit()
    elif nv.EdepMethod =="BetheBloch":
        # Bethe-Bloch dE/dx converted to energy deposit
        # factor 100 is needed to convert wWidth in [m] to [cm]
        nv.enemat = bb.BetheBloch(nv.Particle,float(nv.BEnergy),nv.Material)*nv.Material.rho*nv.WIRESCAN_wWidth*100*np.pi/4 
    else:
        print("Unknown DEDx method: verify!")    

    
    # 
    # If there are electrons in the incident particle the code will also look for the 
    # electron energy deposition file. Be aware that the program will directly look for this file 
    # in the EneDepData folder, without asking. 
    #

    if nv.Particle.Nelectrons != 0:
        Energy_Particle_Electron = float(nv.BEnergy)*nv.Emass/nv.Particle.PartMass

        EleEneDep_filename = "EneDepData/"+"Electron_"+nv.Material.name+".txt"

        h = open(EleEneDep_filename)
      
        v_EleEneDep = []; cont = 0
        for l in h:
            if cont == 0:
                cont += 1; continue
            else:
                v_EleEneDep += [float(l.split()[1])]
                
                if Energy_Particle_Electron <= v_EleEneDep[-1]:
                    if len(v_EleEneDep) <= 1:
                        nv.Ele_enemat = v_EleEneDep[-1]
                        break
                    else: 
                 
                        nv.Ele_enemat = (v_EleEneDep[-2]+v_EleEneDep[-1])/2.0
                        break

        if nv.Ele_enemat == 0:
            nv.Ele_enemat = v_EleEneDep[-1]
        
        h.close()

    #    
    # 6) Load Beam distribution if nongaussian beam. 
    # 
    # If a non gaussian beam type was selected, we will need an input file 
    # indicating how many particles reach each point in space. 
    # See the BeamDescription folder for examples of this input files. 
    # 

    if nv.BeamType != "Gaussian":
        BeamShape_Filename = nv.BeamType
        p = open(BeamShape_Filename)
        nv.Mat_BeamShape = []
        for i,l in enumerate(p,0):
            if i > 1:
                nv.Mat_BeamShape.append([float(l.split()[0]),float(l.split()[1]),float(l.split()[2])])
        p.close()


        
        
        
        
# -------------     Write Output File    ------------------------------------------------------------- #
#
# This function is executed after the simulation is finished, and it creates a series of output 
# files with the requested information. Notice that if one wants other information to be stored 
# this can be done by modifiying this function and adding the corresponding lines. 

def WriteOutputPlotsTxt(foldername):
    

    #  ----- Write Maximum temperature vs time. ------ #

    f1 = open(foldername+"MaxTempVSTime.txt","w")
    f1.write("# ------------ Maximum Temperature Vs Time ----------- #\n")
    f1.write("# input file: "+nv.RealInputFilename+" \n")
    if (nv.DetType == "SEM") or (nv.DetType == "FOIL") or (nv.DetType == "SPLITTER"):
        f1.write("#   Emissivity   |   Time [us]  |   Temperature  [K] |   \n")
    else: 
        f1.write("#   Emissivity  |  Time [us]  |  Position [cm]  |  Temperature [K] | \n")
     
    for j in range(0,len(nv.V_MaximumTemperature)):
        if (nv.DetType == "SEM") or (nv.DetType == "FOIL") or (nv.DetType == "SPLITTER"):
            f1.write(str(nv.V_Emissivity[j])+"   "+str(round(nv.V_Time[j]*1e+6,6))+"   "+str(round(nv.V_MaximumTemperature[j],3))+"\n")
        else: 
            f1.write(str(nv.V_Emissivity[j])+"   "+str(round(nv.V_Time[j]*1e+6,6))+"   "+str(round(nv.V_Pos[j]*1e+2,6))+"   "+str(round(nv.V_MaximumTemperature[j],3))+"\n")

    f1.close()
        
    #  ----- Write Maximum current vs time. ------ #

    f2 = open(foldername+"IntensityVSTime.txt","w")
    f2.write("# ------------ Intensity Vs Time ----------- #\n\n")
    if (nv.DetType == "SEM") or (nv.DetType == "SPLITTER") : 
        f2.write("#   Time [us]  |  WirePos [cm]   | Intensity [mA] #\n")

    elif nv.DetType == "FOIL":
        f2.write("#   Time [us]  |  Intensity [mA] #\n")

    else: 
        f2.write("#   Time [us]  |   Position [cm]    |   Intensity  [mA] #\n")
      
    for j in range(0,len(nv.V_Time)):
        if nv.DetType == "SEM" or nv.DetType == "SPLITTER":
            f2.write(str(nv.V_Time[j]*1e+6)+"   [ ")
            for nw in range(0,len(nv.xvec)):
                f2.write(str(nv.xvec[nw])+", ")
            f2.write(" ]    [ ")
            for n in range(0,len(nv.V_Current2)):
                f2.write(str(nv.V_Current2[n][j]*1e3)+", ")
            f2.write(" ] \n")
        elif nv.DetType == "FOIL":
            f2.write(str(nv.V_Time[j]*1e+6)+"   "+str(nv.V_Current2[j])+"\n")
        else: 
            f2.write(str(nv.V_Time[j]*1e+6)+"   "+str(round(nv.V_Pos[j]*1e+2,3))+"    ")
            f2.write(str(nv.V_Current2[j]*1e3)+"   ")
            f2.write("\n")

    f2.close()



# -------------     Write Output File    ------------------------------------------------------------- #
#
# This is a new function which outputs the data, 2024.03.29 mariusz.sapinski@psi.ch
# in future user should be able to choose what he wants to save
# and the header should contain main simulation parameters (beam size, intensity etc)

def WriteResults(foldername):
   
    f1 = open(foldername+"Output.txt","w")
    f1.write("# --------------- PyTT output file ------------------ \n")
    f1.write("# input file: "+nv.RealInputFilename+" \n")
    f1.write("# execution date and time: "+datetime.datetime.now().strftime('%Y/%m/%d %H:%M:%S')+" \n")
    f1.write("# --------------------------------------------------- \n")
    if (nv.DetType == "SEM") or (nv.DetType == "FOIL") or (nv.DetType == "SPLITTER"):
        f1.write("#   Emissivity   |   Time [us]  |   Temperature  [K] |   \n")   # correct
    else:   # WIRESCANNER
        f1.write("#   Time [us], Position [mm], Npart, MaxTemp [K], SEMcurr [uA], TOTcurr [uA] \n")
        
    for j in range(0,len(nv.V_MaximumTemperature)):
        if (nv.DetType == "SEM") or (nv.DetType == "FOIL") or (nv.DetType == "SPLITTER"):
          f1.write(str(nv.V_Emissivity[j])+"   "+str(round(nv.V_Time[j]*1e+6,6))+"   "+str(round(nv.V_MaximumTemperature[j],3))+"\n")
        else:  # WIRESCANNER
          # formatting output:
          otime=str(round(nv.V_Time[j]*1e+6,6))   # time converted from s to us
          opos=+str(round(nv.V_Pos[j]*1e+3,6))    # position converted from m to mm
          onpa=+str()
          f1.write(otime+","+opos+","+str(round(nv.V_MaximumTemperature[j],3))+"\n")

    f1.close()        
