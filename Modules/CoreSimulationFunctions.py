# --------------------- Core Simulation functions ---------------------- # 
# 
#  In this file is where all the magic happens. The different functions execute the simulations 
#  for the different detectors. 
#

import numpy as np
import sys, os
from Modules import TargetGeometry
from Modules import NecessaryVariables as nv
from Modules import TempPhysicalModels


# ------------  SEM GRID SIMULATION ------------------------ #
# in PSI it is called harfe monitor
# 

def TempEvolSEM():
    
    # In the input file the user can give, either the number of particles or the Intensity and beam pulse length. 
    # Here we calculate the Number of particles if it was not given. 
    
    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*nv.tpulse/nv.Qe
        
        
    # The time has been discretized. Because the beam time can be very small compared with the cooling time, 
    # We will have two different time steps. One for the "beam" and one for teh cooling time. 
    #
    # We divide the simulation in beam shots. Each beam shot will have a beam time and a cooling time. 
    # the total simulation time will depend on how many beam shots there are and how many 
    # time steps we take per shot. 
    
    numberStepPulse = round(nv.tpulse/nv.dtPulse)
    numberStepCooling = round((1./nv.frec - nv.tpulse)/nv.dtCooling)


    if (numberStepPulse == 0) or (numberStepPulse == 0):
        print("dtPulse or dtCooling have been uncorrectly choosen! ")
        sys.exit()

    totalSteps = numberStepCooling + numberStepPulse
    TOTAL = totalSteps*nv.Npulses
    
    # 
    # Here we initialize the necessary variables that will be updated during the simulation. 
    # 

    Time = np.ones(TOTAL); t = 0; Time[0] = t                 # Time of the simulation step
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])       # Matrix with the temperature distribution at each time step. 
    #Tmax = np.ones(TOTAL); Tmax[0] = Temp.max()               # Maximum temperature at each time step. 
    Tmax = np.ones([nv.SEM_nWires,TOTAL]); #Tmax[0] = Temp.max()               # Maximum temperature at each time step. 
    
    
    # Notice that saving all the time the temperature and intensity distribution might slow down the code and increase significantly the output files. 
    FancyI = 0.0*np.ones([nv.SEM_nWires,TOTAL])               # This records the information of the Intensity distribution for all the time steps. 
    FancyT = 0.0*np.ones([nv.SEM_nWires,TOTAL])               # This records the information of the temperature distribution for all the time steps. 
  
    Imax1 = 0.0*np.ones([nv.SEM_nWires,TOTAL])
    Imax2 = 0.0*np.ones([nv.SEM_nWires,TOTAL])
    #Imax1 = 0.0*np.ones([TOTAL])
    #Imax2 = 0.0*np.ones([TOTAL])



    Flag_Current = 0.0           #   This is just to calculate or not the current and store it. 

    l = 0       #   Step Counter

    #    A cycle indicates the beam pulse number. 
    #    A step indicates the time step in that beam pulse. 

    for cycle in range(0,nv.Npulses):
        for step in range(0,int(totalSteps)):
            print("pulse = ",cycle," step = ",step," out of ",totalSteps)

            #  If the material properties are set to vary with the temperature they need to be updated at each 
            #  time step as the temperature will change. 
            #  If the material properties are set to be constant we consider in each time step a temperature 
            #  of 300 K whe checking the value of the parameter. 
            #  
            #  Notice that this gives a matrix for each parameter value, with the same size as the detector geometry.
            #  For each detector position, there is 
            #  a different temperature and thus a different parameter value. 
            
            if nv.EnableParameterVariation:
                nv.Material.epsT = nv.Material.GetEmissivity(Temp)
                nv.Material.CpT = nv.Material.GetCp(Temp)
                nv.Material.wfun = nv.Material.GetWf(Temp)                
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
                nv.Material.wfun = nv.Material.GetWf(300*Temp**0)                
                nv.Material.con = nv.Material.Getk(300*Temp**0)
                nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
            # this stores the value of the emissivity at each time step. It can be removed if it is not of interest. 
            # 2do: either remove or do the same for other temperature-dependent parameters
            nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]
            
            #
            #  For the first time step of the first beam shot we only consider heating as the detector should be at the 
            # given initial temperature. 
            #
            if cycle == 0 and step == 0:
                cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0
                #heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                #heat,_ = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                dt=nv.dtPulse
                #heat,hene,npart = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                heat,hene,npart = TempPhysicalModels.BeamHeating(Temp, dt)
                
                nv.V_Npar.append(npart)


                print("Temp=",Temp)
                print("heat=",heat)
                Temp = Temp + heat

                Flag_Current = 1.0
                
                t += nv.dtPulse

                # This is also used for calculating useful output parameters. We find the index of the highest point in the 
                # temperature matrix in order to afterwards obtain the cooling at that point, the different cooling processes 
                # contributions at that point etc. 
                # It can be deleted if it is not of interest. 
                iMaxTemp = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)

                
            #   
            # Here we are at time steps during the beam pulse. Meaning there is a heatng contribution. 
            #
            elif step <= numberStepPulse:
                
                #heat, hene, npart = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                heat, hene, npart = TempPhysicalModels.BeamHeating(Temp, dt)                
                print("numberStepPulse = ", numberStepPulse) # Debug
                cold1,cene1 = TempPhysicalModels.RadiativeCooling(nv.dtPulse, Temp)
                cold2,cene2 = TempPhysicalModels.ThermionicCooling(nv.dtPulse,Temp)
                
                print("npart=",npart)
                
                nv.V_Npar.append(npart)

                
                # Because conduction and sublimation cooling require quite intense calculations, they are not calculated 
                # unless it has been specified. 
                
                if nv.ConductiveCooling == 1:
                    cold3, cene3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else: cold3 = Temp*0.0; cene3 = Temp*0.0

                if nv.SublimationCooling == 1:
                    cold4,cene4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0; cene4=Temp*0.0

                # Depending on the users choices, some cooling effects might not be considered.
                # this nv.RadiativeCooling, nv.ThermioniCooling, etc. Are flags (0: Not Active), (1: Active)
                #print("Temp=",Temp)
                #print("heat=",heat)
                #Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4
                # or differently, add up energies and calculate temperature after
                # 2025.01.26
                dEne=hene-cene1-cene2-cene3-cene4
                Temp = Temp + dEne/(nv.Material.CpT*nv.eVol*nv.Material.rho*1e+6)

                Flag_Current = 1.0

                t += nv.dtPulse
                
            # 
            # Here we are at the point where there is no more beam pulse anymore, so the heating term is not there. 
            # Also here the time step is the one for the cooling part. 
            # 
            else:
                cold1,cene1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2,cene2 = TempPhysicalModels.ThermionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3,cene3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0; cene3 = Temp * 0.0 

                if nv.SublimationCooling == 1:
                    cold4,cene4 = TempPhysicalModels.SublimationCooling(nv.dtCooling, Temp)
                else: cold4 = Temp * 0.0; cene4 = Temp * 0.0

                #Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4
                # 2025.01.26
                dEne=hene-cene1-cene2-cene3-cene4
                Temp = Temp + dEne/(nv.Material.CpT*nv.eVol*nv.Material.rho*1e+6)



                Flag_Current = 0.0      # The current is not clculated in this part of the simulation. Only when there is beam. It can be changed. 

                t += nv.dtCooling

                totalcooling = cold1[iMaxTemp]*nv.RadiativeCooling + cold2[iMaxTemp]*nv.ThermionicCooling + cold3[iMaxTemp]*nv.ConductiveCooling+nv.SublimationCooling*cold4[iMaxTemp]
                nv.CoolingImportance_Temp += [Temp[iMaxTemp]]
                
                nv.CoolingImportance_Ems += [cold1[iMaxTemp]*nv.RadiativeCooling/totalcooling]
                nv.CoolingImportance_Jth += [cold2[iMaxTemp]*nv.ThermionicCooling/totalcooling]
                nv.CoolingImportance_Con += [cold3[iMaxTemp]*nv.ConductiveCooling/totalcooling]
                nv.CoolingImportance_Sub += [cold4[iMaxTemp]*nv.SublimationCooling/totalcooling]
            
            l += 1
            
            
            #
            #  ------------ Calculate Current Generated in Wire ------------ #
            # Notice that by default this current is only calculated during the beam pulse. 
            # 
            # before 2025.01.29:
            #V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
            #current2 = V_current[1]
            #  ------------ Calculate Current Generated in Wire ------------ #
            V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
            current1 = V_current[0]  # Secondary electron current [A]   # it contains all wires!
            current2 = V_current[1]  # Thermionic emission current [A]
            # -------------------------------------------------------------- #            
            print(type(l),l, current1)
            print(type(current1))
            print(type(Imax1))
            print(Imax1.shape)
            for iwir in range(nv.SEM_nWires): 
                Imax1[iwir][l-1] = current1[iwir]   # Imax1 is a misleading name I guess
                Imax2[iwir][l-1] = current2[iwir]            


            #
            # In this part we store all the useful information that was colected in this time step. 

            if l > nv.Npulses * totalSteps - 1:
                continue
            else:
                Time[l] = t
                #print("Debug, Temperature = ",Temp)
                for iwir in range(nv.SEM_nWires):
                    #print("Debug, Temperature = ",iwir,Temp[iwir],np.max(Temp[iwir]))
                    Tmax[iwir][l-1] = np.max(Temp[iwir][:])
                    #print(Tmax[iwir][l-1])

                for mk in range(0,len(current2)):
                    FancyI[mk][l] = current2[mk]

                # We need to check if it is a horizontal or vertical grid. 
                for mt in range(0,nv.SEM_nWires):
                    if nv.SEM_Plane == "Horizontal":
                        FancyT[mt][l] = np.max(Temp[mt][:])
                    else: 
                        FancyT[mt][l] = np.max(Temp[:][mt])
                        
                        
            print("Simulation: ", l, "    From: ", nv.Npulses*totalSteps, "    Tmax: ", np.max(Temp))

 


            
    #return Time, Tmax, FancyTemp, Imax1, V_Pos, Imax2  # from WireSCAN1
    #        0     1      2      3     
    # old (before 2025.01.29):
    # return Time, Tmax, FancyT, FancyI
    # Imax1, Imax2 - arrays of currents, for SEM - per wire
    return Time, Tmax, FancyT, Imax1, Imax2


# ------------  SPLITTER SIMULATION ------------------------ #
# how this is different from SEM GRID? Only wire shape?

def TempEvolSPLITTER():
    
    # In the input file the user can give, either the number of particles or the Intensity and beam pulse length. 
    # Here we calculate the Number of particles if it was not given. 
    
    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*nv.tpulse/nv.Qe
        
        
    # The time has been discretized. Because the beam time can be very small compared with the cooling time, 
    # We will have two different time steps. One for the "beam" and one for teh cooling time. 
    #
    # We divide the simulation in beam shots. Each beam shot will have a beam time and a cooling time. 
    # the total simulation time will depend on how many beam shots there are and how many 
    # time steps we take per shot. 
    
    numberStepPulse = round(nv.tpulse/nv.dtPulse)
    numberStepCooling = round((1./nv.frec - nv.tpulse)/nv.dtCooling)  # what is that?

    if (numberStepPulse == 0) or (numberStepCooling == 0):
        print("dtPulse or dtCooling have been uncorrectly choosen! ")
        sys.exit()

    totalSteps = numberStepCooling + numberStepPulse
    TOTAL = totalSteps
    
    # 
    # Here we initialize the necessary variables that will be updated during the simulation. 
    # 

    Time = np.ones(TOTAL); t = 0; Time[0] = t                 # Time of the simulation step
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])       # Matrix with the temperature distribution at each time step. 
    Tmax = np.ones(TOTAL); Tmax[0] = Temp.max()               # Maximum temperature at each time step. 
    
    # Notice that saving all the time the temperature and intensity distribution might slow down the code and increase significantly the output files. 
    #FancyI = 0.0*np.ones([nv.SEM_nWires,TOTAL])               # This records the information of the Intensity distribution for all the time steps. 
    #FancyT = 0.0*np.ones([nv.SEM_nWires,TOTAL])               # This records the information of the temperature distribution for all the time steps. 
  

    Flag_Current = 0.0           #   This is just to calculate or not the current and store it, not used

    l = 0       #   Step Counter

    #    A cycle indicates the beam pulse number. 
    #    A step indicates the time step in that beam pulse. 

    #for cycle in range(0,nv.Npulses):
    print("total steps = ",totalSteps)
    for step in range(0,int(totalSteps)):

        print("next step nr ",step)
        #  If the material properties are set to vary with the temperature they need to be updated at each 
        #  time step as the temperature will change. 
        #  If the material properties are set to be constant we consider in each time step a temperature 
        #  of 300 K whe checking the value of the parameter. 
        #  
        #  Notice that this gives a matrix for each parameter value, with the same size as the detector geometry.
        #  For each detector position, there is 
        #  a different temperature and thus a different parameter value. 
            
        if nv.EnableParameterVariation:
            nv.Material.epsT = nv.Material.GetEmissivity(Temp)
            nv.Material.CpT = nv.Material.GetCp(Temp)
            nv.Material.wfun = nv.Material.GetWf(Temp)            
            nv.Material.con = nv.Material.Getk(Temp)
            nv.Material.HT = nv.Material.GetH(Temp)
        else: 
            nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
            nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
            nv.Material.wfun = nv.Material.GetWf(300*Temp**0)
            nv.Material.con = nv.Material.Getk(300*Temp**0)
            nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
        # this stores the value of the emissivity at each time step. It can be removed if it is not of interest. 
        nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]
            
        #
        #  For the first time step of the first beam shot we only consider heating as the detector should be at the 
        # given initial temperature. 
        #
        if step == 0: # first step
            cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0
            # call this once to get SEYp info
            Npart=np.ones(TOTAL)
            dt=1e-3
            Current1, Current2 = TempPhysicalModels.CalculateCurrent(Npart,Temp,numberStepPulse,dt)
            heat,_ = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
            Temp = Temp + heat

            Flag_Current = 1.0
                
            t += nv.dtPulse

            # This is also used for calculating useful output parameters. We find the index of the highest point in the 
            # temperature matrix in order to afterwards obtain the cooling at that point, the different cooling processes 
            # contributions at that point etc. 
            # It can be deleted if it is not of interest. 
            iMaxTemp = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)

                
        #   
        # Here we are at time steps during the beam pulse. Meaning there is a heatng contribution. 
        #
        elif step <= numberStepPulse:
                
            heat,_ = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
            cold1 = TempPhysicalModels.RadiativeCooling(nv.dtPulse, Temp)
            cold2 = TempPhysicalModels.ThermionicCooling(nv.dtPulse,Temp)
                
            # Because conduction and sublimation cooling require quite intense calculations, they are not calculated 
            # unless it has been specified. 
                
            if nv.ConductiveCooling == 1:
                cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
            else: cold3 = Temp*0.0


            # Depending on the users choices, some cooling effects might not be considered.
            # this nv.RadiativeCooling, nv.ThermioniCooling, etc. Are flags (0: Not Active), (1: Active)
            Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3

            Flag_Current = 1.0

            t += nv.dtPulse
                
            # 
            # Here we are at the point where there is no more beam pulse anymore, so the heating term is not there. 
            # Also here the time step is the one for the cooling part. 
            # 
        else:
            cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
            cold2 = TempPhysicalModels.ThermionicCooling(nv.dtCooling,Temp)
            if nv.ConductiveCooling == 1:
                cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
            else:
                cold3 = Temp * 0.0

            Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3

            Flag_Current = 0.0      # The current is not clculated in this part of the simulation. Only when there is beam. It can be changed. 

            t += nv.dtCooling

            totalcooling = cold1[iMaxTemp]*nv.RadiativeCooling + cold2[iMaxTemp]*nv.ThermionicCooling + cold3[iMaxTemp]*nv.ConductiveCooling
            nv.CoolingImportance_Temp += [Temp[iMaxTemp]]
                
            nv.CoolingImportance_Ems += [cold1[iMaxTemp]*nv.RadiativeCooling/totalcooling]
            nv.CoolingImportance_Jth += [cold2[iMaxTemp]*nv.ThermionicCooling/totalcooling]
            nv.CoolingImportance_Con += [cold3[iMaxTemp]*nv.ConductiveCooling/totalcooling]
            
        l += 1  # what is l?
            
            
        #
        #  ------------ Calculate Current Generated in Wire ------------ #
        # Notice that by default this current is only calculated during the beam pulse. 
        # 
        #V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
        #current2 = V_current[1]
            
           
        #
        # In this part we store all the useful information that was colected in this time step. 

        #if l > nv.Npulses * totalSteps - 1:
        #    continue
        #else:
        #    Time[l] = t
        #    Tmax[l] = np.max(Temp)
        Time[step] = t
        Tmax[step] = np.max(Temp)

                        
                        
        print("Simulation: ", step, "    From: ", totalSteps, "    Tmax: ", np.max(Temp))
            
            
    return Time, Tmax





# ------------------------------------------ FOIL ---------------------------------------------- #
# 
# The simulations for the foil are basically identical to the ones for the SEM grid. Please 
# checkk the coments on the TimeEvolSEM() function. 

def TimeEvolFOIL():

    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*nv.tpulse/nv.Qe

    numberStepPulse = round(nv.tpulse/nv.dtPulse)
    numberStepCooling = round((1./nv.frec - nv.tpulse)/nv.dtCooling)

    if (numberStepPulse == 0) or (numberStepPulse == 0):
        print("dtPulse or dtCooling have been uncorrectly choosen! ")
        sys.exit()


    totalSteps = numberStepCooling + numberStepPulse
    TOTAL = totalSteps*nv.Npulses

    Time = np.ones(TOTAL); t = 0; Time[0] = t
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])
    Tmax = np.ones(TOTAL); Tmax[0] = Temp.max()

    Imax1 = 0.0*np.ones(TOTAL)
    Imax2 = 0.0*np.ones(TOTAL)
    Flag_Current = 0.0

    l = 0       # Counter

    for cycle in range(0,nv.Npulses):
        for step in range(0,int(totalSteps)):
            if nv.EnableParameterVariation:
                nv.Material.epsT = nv.Material.GetEmissivity(Temp)
                nv.Material.CpT = nv.Material.GetCp(Temp)
                nv.Material.wfun = nv.Material.GetWf(Temp)                
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
                nv.Material.wfun = nv.Material.GetWf(300*Temp**0)
                nv.Material.con = nv.Material.Getk(300*Temp**0)
                nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
            nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]

            if cycle == 0 and step == 0:
                cold1 = 0; cold2 = 0; cold3 = 0
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                Temp = Temp + heat

                Flag_Current = 1.0
                
                t += nv.dtPulse

                iMaxTemp = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)

            elif step <= numberStepPulse:
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtPulse, Temp)
                cold2 = TempPhysicalModels.ThermionicCooling(nv.dtPulse,Temp)
                
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else:
                    cold3 = Temp*0.0

                Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3

                Flag_Current = 1.0

                t += nv.dtPulse

            else:
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2 = TempPhysicalModels.ThermionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0
                Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3

                Flag_Current = 0.0

                t += nv.dtCooling

                totalcooling = cold1[iMaxTemp]*nv.RadiativeCooling + cold2[iMaxTemp]*nv.ThermionicCooling + cold3[iMaxTemp]*nv.ConductiveCooling
                nv.CoolingImportance_Temp += [Temp[iMaxTemp]]
                
                nv.CoolingImportance_Ems += [cold1[iMaxTemp]*nv.RadiativeCooling/totalcooling]
                nv.CoolingImportance_Jth += [cold2[iMaxTemp]*nv.ThermionicCooling/totalcooling]
                nv.CoolingImportance_Con += [cold3[iMaxTemp]*nv.ConductiveCooling/totalcooling]
                
                
            
            l += 1
            #if l == 18: sys.exit()
            print("Simulation: ", l, "    From: ", nv.Npulses*totalSteps, "    Tmax: ", np.max(Temp))
            
            #  ------------ Calculate Current Generated in Wire ------------ #
            V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
            current1 = V_current[0]
            current2 = V_current[1]

            if l > nv.Npulses * totalSteps - 1:
                continue
            else:
                Time[l] = t
                Tmax[l] = np.max(Temp)
                Imax1[l] = current1
                Imax2[l] = current2
                

            # This is for temperature distribution plot. 
            if l == nv.Npulses*totalSteps - numberStepCooling:
                FancyTemp = Temp

    nv.WireExp = TempPhysicalModels.LinearThermalExpansion(FancyTemp)

    return Time, Tmax, FancyTemp, Imax1, Imax2


# --------------------------------- SLOW WIRE SCANNER -------------------------------------- #
# 
# In the case of a slow wire scanner, the scanner position will be updated for every beam pulse. 
# That is, at the beggining of a beam pulse the position of the scaner will be moved the amount of distance
# that was specified in the input. 
# aside from this update in the position the simulation is very similar to the case of the 
# foils and the SEM grids. 
#

def TimeEvolWIRESCAN2():

    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*nv.tpulse/nv.Qe

    numberStepPulse = round(nv.tpulse/nv.dtPulse)
    numberStepCooling = round((1./nv.frec - nv.tpulse)/nv.dtCooling)

    if (numberStepPulse == 0) or (numberStepPulse == 0):
        print("dtPulse or dtCooling have been uncorrectly choosen! ")
        sys.exit()

    totalSteps = numberStepCooling + numberStepPulse
    TOTAL = totalSteps*nv.Npulses

    # 
    #  Here we need to calculate first which postions the wire will be taking. The pos_vector stores this 
    # positions. They will be calculated depending on the initial/final position specified and nv.Npulses
    # will determine how many steps there are in between this initial, final position. 
    # Note: This could be improved, we could just give a number of points but for the moment it is like this. 
    # 
    
    pos_vector = []
    
    
    if nv.WIRESCAN_Plane == "Horizontal":
        x0 = nv.WIRESCAN_IniPos; x1 = nv.WIRESCAN_EndPos
        for k in range(0,nv.Npulses+1):
            pos_vector += [x0+(x1-x0)*k/nv.Npulses]
        nv.WIRESCAN_wCposX = pos_vector[0]
        nv.WIRESCAN_wCposY = 0.0
        TargetGeometry.CreateDetector(nv.DetType)

    elif nv.WIRESCAN_Plane == "Vertical":
        y0 = nv.WIRESCAN_IniPos; y1 = nv.WIRESCAN_EndPos
        for k in range(0,nv.Npulses+1):
            pos_vector += [y0+(y1-y0)*k/nv.Npulses]
        nv.WIRESCAN_wCposX = 0.0
        nv.WIRESCAN_wCposY = pos_vector[0]
        TargetGeometry.CreateDetector(nv.DetType)

    
    
    # 
    # As before we initialize the necessary variables. 
    # 

    Time = np.ones(TOTAL); t = 0; Time[0] = t
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])
    Tmax = np.ones(TOTAL); Tmax[0] = Temp.max()

    Imax1 = 0.0*np.ones(TOTAL)
    Imax2 = 0.0*np.ones(TOTAL)

    V_Pos = 0.0*np.ones(TOTAL)

    Flag_Current = 0.0

    l = 0       # Counter
    for cycle in range(0,nv.Npulses):
        
        # ------ Update Wire Position and Particle Matrix ------ #
        # We need to update the positon but also how many particles reach the wire at that 
        # new position. Here the wire "geometry" is created at each wire step.
        #

        if nv.WIRESCAN_Plane == "Horizontal":
            nv.WIRESCAN_wCposX = pos_vector[cycle]
            nv.WIRESCAN_wCposY = 0.0
            TargetGeometry.CreateDetector(nv.DetType)
        elif nv.WIRESCAN_Plane == "Vertical":
            nv.WIRESCAN_wCposX = 0.0
            nv.WIRESCAN_wCposY = pos_vector[cycle]
            TargetGeometry.CreateDetector(nv.DetType)
        elif nv.WIRESCAN_Plane == "Diagonal":
            print("Diagonally Moving wire scanner not ready yet sorry!")
            sys.exit()

        nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()
        heat = 0; cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0
        # ------------------------------------- #

        for step in range(0,int(totalSteps)):
            if nv.EnableParameterVariation:
                nv.Material.epsT = nv.Material.GetEmissivity(Temp)
                nv.Material.CpT = nv.Material.GetCp(Temp)
                nv.Material.wfun = nv.Material.GetWf(Temp)                
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
                nv.Material.wfun = nv.Material.GetWf(300*Temp**0)                
                nv.Material.con = nv.Material.Getk(300*Temp**0)
                nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
            nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]
            if cycle == 0 and step == 0:
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                Temp = Temp + heat
                Flag_Current = 1.0
                t += nv.dtPulse

            elif step <= numberStepPulse:
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtPulse, Temp)
                cold2 = TempPhysicalModels.ThermionicCooling(nv.dtPulse,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else:
                    cold3 = Temp*0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0

                
                Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3+nv.SublimationCooling*cold4

                Flag_Current = 1.0

                t += nv.dtPulse
            else:
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2 = TempPhysicalModels.ThermionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0

                Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3+nv.SublimationCooling*cold4

                Flag_Current = 0.0

                t += nv.dtCooling

                nv.CoolingImportance_Temp += [np.max(Temp)]
                nv.CoolingImportance_Ems += [np.min(cold1*nv.RadiativeCooling)]
                nv.CoolingImportance_Jth += [np.min(cold2*nv.ThermionicCooling)]
                nv.CoolingImportance_Con += [np.min(cold3*nv.ConductiveCooling)]
           
            l += 1; 
            print("Simulation: ", l, "    From: ", nv.Npulses*totalSteps, "    Tmax: ", np.max(Temp), "Cycle: ",cycle, "Pos: ", nv.xvec, " Part: ",np.max(nv.ParticleProportionMatrix)*nv.Nparticles)
            #  ------------ Calculate Current Generated in Wire ------------ #
            V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
            current1 = V_current[0]
            current2 = V_current[1]            
            # -------------------------------------------------------------- #
            if l > nv.Npulses * totalSteps - 1:
                continue
            else:
                Time[l] = t
                Tmax[l] = np.max(Temp)
                
                Imax1[l] = np.max(current1)  # should it not be sum?
                Imax2[l] = np.max(current2)

                if nv.WIRESCAN_Plane == "Horizontal": 
                    V_Pos[l] = nv.xvec[0] 
                elif nv.WIRESCAN_Plane == "Vertical": 
                    V_Pos[l] = nv.yvec[0]

            # This is for temperature distribution plot. 
            if l == nv.Npulses*totalSteps - numberStepCooling:
                if nv.WIRESCAN_Plane == "Horizontal":
                    FancyTemp = Temp
                elif nv.WIRESCAN_Plane == "Vertical":
                    FancyTemp = np.transpose(Temp)

    nv.WireExp = TempPhysicalModels.LinearThermalExpansion(FancyTemp)

    return Time, Tmax, FancyTemp, [Imax2], V_Pos



# ---------------------------- FAST WIRE SCANNER ------------------------------------- #
# 
#  This simulation is slightly different than the previous ones. It is difficult to match the real 
# beam pulse structure and its movement. For that reason we consider an equivalent constant beam and a wire that moves  
# at a certain speed throught this constant beam. 
# 


def TimeEvolWIRESCAN1():
    
    # 
    #   Here the simulation steps are set by default. The user has no say in this unless this hard coded parameter is modified. 
    # 
    Nsteps = 5000             # Simulations accuracy. How many points we will divide the space in. 
    # idea: put it into configuration file    
    
    # remove files used for debugging in order not to keep them growing:
    if os.path.isfile("Output/Nimatrix.txt"):
        os.remove("Output/Nimatrix.txt")    
    
    if os.path.isfile("Output/beam_profile.txt"):
        os.remove("Output/beam_profile.txt")    


    
   
    if nv.WIRESCAN_wSpeed <= 0.0:
        print("Wire Scanner Velocity negative or zero - not Valid!")
        sys.exit()
        
    #
    #  First we create the wire geometry as in the prevous cases. In this case, we use the 
    # initial position/final position parameters, the wire speed and the number of steps of the simulation 
    # to calculate the different positions of the wire that will be simulated. 
    #
    
    if nv.WIRESCAN_Plane == "Vertical":
        y0 = nv.WIRESCAN_IniPos; y1 = nv.WIRESCAN_EndPos
        nv.WIRESCAN_wCposX = 0.0
        nv.WIRESCAN_wCposY = y0
        dt = abs(y1-y0)/(nv.WIRESCAN_wSpeed*Nsteps)

    elif nv.WIRESCAN_Plane == "Horizontal":
        x0 = nv.WIRESCAN_IniPos; x1 = nv.WIRESCAN_EndPos
        nv.WIRESCAN_wCposX = x0
        nv.WIRESCAN_wCposY = 0.0
        dt = abs(x1-x0)/(nv.WIRESCAN_wSpeed*Nsteps)   # [s] calculate recommended dt
        
    if dt>nv.WIRESCAN_wWidth/nv.WIRESCAN_wSpeed:
        print("Problem with too long dt: ",dt)
        print(">>should be at least: ",nv.WIRESCAN.wWidth/nv.WIRESCAN_wSpeed)
        exit(0)
        
        
    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*dt/nv.Qe
 
    # ms: new May 19, 2024, use beam intensity [A] as principal parameter
    if nv.Intensity == 0.0:
        nv.Intensity = nv.Nparticles * nv.frec * nv.Qe  # here protons are assumed 
        
 
    
 
    TargetGeometry.CreateDetector(nv.DetType)
    
    #
    # We initialize the necesary parameters
    # 
    
    Time = np.ones(Nsteps); t = 0; Time[0] = t
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])
    Tmax = np.ones(Nsteps); Tmax[0] = Temp.max()

    Imax1 = 0.0*np.ones(Nsteps)
    Imax2 = 0.0*np.ones(Nsteps)

    V_Pos = 0.0*np.ones(Nsteps)

    Flag_Current = 1.0

    #nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()  # not needed anymore

    heat = 0; cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0

    for step in range(0,Nsteps):
       
        # 
        # Beacause of the movement at each  time steps we calculate the amount of particles reaching the wire. 
        #
        
        nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()  # do we need it?
        nv.stepcount+=1

        if nv.EnableParameterVariation:
            nv.Material.epsT = nv.Material.GetEmissivity(Temp)
            nv.Material.CpT = nv.Material.GetCp(Temp)
            nv.Material.wfun = nv.Material.GetWf(Temp)            
            nv.Material.con = nv.Material.Getk(Temp)
            nv.Material.HT = nv.Material.GetH(Temp)
        
        else: 
            nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
            nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
            nv.Material.wfun = nv.Material.GetWf(300*Temp**0)
            nv.Material.con = nv.Material.Getk(300*Temp**0)
            nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
        nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]

        heat,hene,npart = TempPhysicalModels.BeamHeating(Temp, dt)
        nv.V_Npar.append(npart)
        #print(heat,npart)
        #sys.exit()
        #print(np.max(heat)); sys.exit()
        if nv.RadiativeCooling == 1:
            cold1,cene1 = TempPhysicalModels.RadiativeCooling(dt, Temp)
        else: cold1 = Temp*0.0; cene1 = Temp*0.0

        if nv.ThermionicCooling == 1:
            cold2,cene2 = TempPhysicalModels.ThermionicCooling(dt,Temp)
        else: cold2 = Temp*0.0; cene2 = Temp*0.0

        if nv.ConductiveCooling == 1:
            cold3,cene3 = TempPhysicalModels.ConductiveCooling(dt, Temp)
        else: cold3 = Temp*0.0; cene3 = Temp*0.0
        
        if nv.SublimationCooling == 1:
            cold4 = TempPhysicalModels.SublimationCooling(dt, Temp)
        else: cold4 = Temp*0.0

        # add up all the temperature changes:
        #Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4

        # or differently, add up energies and calculate temperature after
        # 2024.05.20
        dEne=hene-cene1-cene2
        Temp = Temp + dEne/(nv.Material.CpT*nv.eVol*nv.Material.rho*1e+6)


        t += dt

        nv.CoolingImportance_Temp += [np.max(Temp)]
        nv.CoolingImportance_Ems += [np.min(cold1*nv.RadiativeCooling)]
        nv.CoolingImportance_Jth += [np.min(cold2*nv.ThermionicCooling)]
        nv.CoolingImportance_Con += [np.min(cold3*nv.ConductiveCooling)]


        print("TimeEvolWIRESCAN1 step Nr: ", step, "    From: ", Nsteps, "    Tmax: ", np.max(Temp), "   x: ", nv.xvec*1e+3 , " [mm]")
        #  ------------ Calculate Current Generated in Wire ------------ #
        V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,1.0,dt)
        current1 = V_current[0]  # Secondary electron current [A]
        current2 = V_current[1]  # Thermionic emission current [A]
        # -------------------------------------------------------------- #
      
        Time[step] = t
        Tmax[step] = np.max(Temp)

        Imax1[step] = current1
        Imax2[step] = current2

        FancyTemp = np.transpose(Temp) 

        # -------------- Update Position -------------- #
        if nv.WIRESCAN_Plane == "Horizontal":
            nv.WIRESCAN_wCposX += nv.WIRESCAN_wSpeed*dt # wire centre position
            V_Pos[step] = nv.xvec[0]
    
        elif nv.WIRESCAN_Plane == "Vertical":
            nv.WIRESCAN_wCposY += nv.WIRESCAN_wSpeed*dt
            V_Pos[step] = nv.yvec[0]
          
        TargetGeometry.CreateDetector(nv.DetType)

    nv.WireExp = TempPhysicalModels.LinearThermalExpansion(FancyTemp)

    return Time, Tmax, FancyTemp, Imax1, V_Pos, Imax2
