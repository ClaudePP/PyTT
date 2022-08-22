import numpy as np
import sys
from Modules import TargetGeometry
from Modules import NecessaryVariables as nv
from Modules import TempPhysicalModels

def TempEvolSEM():
    
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

    FancyI = 0.0*np.ones([nv.SEM_nWires,TOTAL])
    FancyT = 0.0*np.ones([nv.SEM_nWires,TOTAL])

    Flag_Current = 0.0

    l = 0       # Counter


    for cycle in range(0,nv.Npulses):
        for step in range(0,int(totalSteps)):

            if nv.EnableParameterVariation:
                nv.Material.epsT = nv.Material.GetEmissivity(Temp)
                nv.Material.CpT = nv.Material.GetCp(Temp)
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
                nv.Material.con = nv.Material.Getk(300*Temp**0)
                nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
            nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]

            if cycle == 0 and step == 0:
                cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                Temp = Temp + heat

                Flag_Current = 1.0
                
                t += nv.dtPulse

                iMaxTemp = np.unravel_index(np.argmax(Temp, axis=None), Temp.shape)

            elif step <= numberStepPulse:
                heat = TempPhysicalModels.BeamHeating(Temp, numberStepPulse)
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtPulse, Temp)
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtPulse,Temp)
                
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else: cold3 = Temp*0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0

                Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4

                Flag_Current = 1.0

                t += nv.dtPulse
                
            else:
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtCooling, Temp)
                else: cold4 = Temp*0.0

                Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4

                Flag_Current = 0.0

                t += nv.dtCooling

                totalcooling = cold1[iMaxTemp]*nv.RadiativeCooling + cold2[iMaxTemp]*nv.ThermoionicCooling + cold3[iMaxTemp]*nv.ConductiveCooling+nv.SublimationCooling*cold4[iMaxTemp]
                nv.CoolingImportance_Temp += [Temp[iMaxTemp]]
                
                nv.CoolingImportance_Ems += [cold1[iMaxTemp]*nv.RadiativeCooling/totalcooling]
                nv.CoolingImportance_Jth += [cold2[iMaxTemp]*nv.ThermoionicCooling/totalcooling]
                nv.CoolingImportance_Con += [cold3[iMaxTemp]*nv.ConductiveCooling/totalcooling]
                nv.CoolingImportance_Sub += [cold4[iMaxTemp]*nv.SublimationCooling/totalcooling]
            
            l += 1
            
            
            
            #  ------------ Calculate Current Generated in Wire ------------ #
            V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,numberStepPulse,nv.dtPulse)
            current2 = V_current[1]
            
           
            # -------------------------------------------------------------- #

            if l > nv.Npulses * totalSteps - 1:
                continue
            else:
                Time[l] = t
                Tmax[l] = np.max(Temp)

                for mk in range(0,len(current2)):
                    FancyI[mk][l] = current2[mk]

                # We need to check if it is a horizontal or vertical grid. 
                for mt in range(0,nv.SEM_nWires):
                    if nv.SEM_Plane == "Horizontal":
                        FancyT[mt][l] = np.max(Temp[mt][:])
                    else: 
                        FancyT[mt][l] = np.max(Temp[:][mt])
            print("Simulation: ", l, "    From: ", nv.Npulses*totalSteps, "    Tmax: ", np.max(Temp))
    return Time, Tmax, FancyT, FancyI

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
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
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
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtPulse,Temp)
                
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else:
                    cold3 = Temp*0.0

                Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3

                Flag_Current = 1.0

                t += nv.dtPulse

            else:
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0
                Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3

                Flag_Current = 0.0

                t += nv.dtCooling

                totalcooling = cold1[iMaxTemp]*nv.RadiativeCooling + cold2[iMaxTemp]*nv.ThermoionicCooling + cold3[iMaxTemp]*nv.ConductiveCooling
                nv.CoolingImportance_Temp += [Temp[iMaxTemp]]
                
                nv.CoolingImportance_Ems += [cold1[iMaxTemp]*nv.RadiativeCooling/totalcooling]
                nv.CoolingImportance_Jth += [cold2[iMaxTemp]*nv.ThermoionicCooling/totalcooling]
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

    # --- Measurement Points ---- #
    # In this simulation type there is one measurement per pulse. 
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

    # --------------------------- #

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
                nv.Material.con = nv.Material.Getk(Temp)
                nv.Material.HT = nv.Material.GetH(Temp)
            else: 
                nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
                nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
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
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtPulse,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtPulse, Temp)
                else:
                    cold3 = Temp*0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0

                
                Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3+nv.SublimationCooling*cold4

                Flag_Current = 1.0

                t += nv.dtPulse
            else:
                cold1 = TempPhysicalModels.RadiativeCooling(nv.dtCooling,Temp)
                cold2 = TempPhysicalModels.ThermoionicCooling(nv.dtCooling,Temp)
                if nv.ConductiveCooling == 1:
                    cold3 = TempPhysicalModels.ConductiveCooling(nv.dtCooling,Temp)
                else:
                    cold3 = Temp * 0.0

                if nv.SublimationCooling == 1:
                    cold4 = TempPhysicalModels.SublimationCooling(nv.dtPulse, Temp)
                else: cold4 = Temp*0.0

                Temp = Temp + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3+nv.SublimationCooling*cold4

                Flag_Current = 0.0

                t += nv.dtCooling

                nv.CoolingImportance_Temp += [np.max(Temp)]
                nv.CoolingImportance_Ems += [np.min(cold1*nv.RadiativeCooling)]
                nv.CoolingImportance_Jth += [np.min(cold2*nv.ThermoionicCooling)]
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
                
                Imax1[l] = np.max(current1)
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

def TimeEvolWIRESCAN1():
    # ---------- Fixed simulation Parameters ------------ #
    Nsteps = 100000             # Simulations accuracy. How many points we will divide the space in. 
    # --------------------------------------------------- #
    if nv.WIRESCAN_wSpeed <= 0.0:
        print("Wire Scanner Velocity not Valid!")
        sys.exit()
    # ---------- Creating Geometry ------------ #
    if nv.WIRESCAN_Plane == "Vertical":
        y0 = nv.WIRESCAN_IniPos; y1 = nv.WIRESCAN_EndPos
        nv.WIRESCAN_wCposX = 0.0
        nv.WIRESCAN_wCposY = y0
        dt = abs(y1-y0)/(nv.WIRESCAN_wSpeed*Nsteps)

    elif nv.WIRESCAN_Plane == "Horizontal":
        x0 = nv.WIRESCAN_IniPos; x1 = nv.WIRESCAN_EndPos
        nv.WIRESCAN_wCposX = x0
        nv.WIRESCAN_wCposY = 0.0
        dt = abs(x1-x0)/(nv.WIRESCAN_wSpeed*Nsteps)

    if nv.Nparticles == 0.0:
        nv.Nparticles = nv.Intensity*dt/nv.Qe
    
    # ------------------------------------------- #
    TargetGeometry.CreateDetector(nv.DetType)
    # ---------- Simulation Parameters ------------ #
    Time = np.ones(Nsteps); t = 0; Time[0] = t
    Temp = nv.T0 * np.ones([len(nv.xvec),len(nv.yvec)])
    Tmax = np.ones(Nsteps); Tmax[0] = Temp.max()

    Imax1 = 0.0*np.ones(Nsteps)
    Imax2 = 0.0*np.ones(Nsteps)

    V_Pos = 0.0*np.ones(Nsteps)

    Flag_Current = 1.0

    nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()

    heat = 0; cold1 = 0; cold2 = 0; cold3 = 0; cold4 = 0

    for step in range(0,Nsteps):
        nv.ParticleProportionMatrix = TempPhysicalModels.CreateNiMatrix()

        if nv.EnableParameterVariation:
            nv.Material.epsT = nv.Material.GetEmissivity(Temp)
            nv.Material.CpT = nv.Material.GetCp(Temp)
            nv.Material.con = nv.Material.Getk(Temp)
            nv.Material.HT = nv.Material.GetH(Temp)
        
        else: 
            nv.Material.epsT = nv.Material.GetEmissivity(300*Temp**0)
            nv.Material.CpT = nv.Material.GetCp(300*Temp**0)
            nv.Material.con = nv.Material.Getk(300*Temp**0)
            nv.Material.HT = nv.Material.GetH(300*Temp**0)
            
        nv.V_Emissivity += [nv.Material.epsT[np.unravel_index(np.argmax(Temp),Temp.shape)]]

        heat = TempPhysicalModels.BeamHeating(Temp, dt)
        #print(np.max(heat)); sys.exit()
        cold1 = TempPhysicalModels.RadiativeCooling(dt, Temp)
        cold2 = TempPhysicalModels.ThermoionicCooling(dt,Temp)

        if nv.ConductiveCooling == 1:
            cold3 = TempPhysicalModels.ConductiveCooling(dt, Temp)
        else: cold3 = Temp*0.0
        
        if nv.SublimationCooling == 1:
            cold4 = TempPhysicalModels.SublimationCooling(dt, Temp)
        else: cold4 = Temp*0.0

        Temp = Temp + heat + nv.RadiativeCooling*cold1 + nv.ThermoionicCooling*cold2 + nv.ConductiveCooling*cold3 + nv.SublimationCooling*cold4

        t += dt

        nv.CoolingImportance_Temp += [np.max(Temp)]
        nv.CoolingImportance_Ems += [np.min(cold1*nv.RadiativeCooling)]
        nv.CoolingImportance_Jth += [np.min(cold2*nv.ThermoionicCooling)]
        nv.CoolingImportance_Con += [np.min(cold3*nv.ConductiveCooling)]


        print("Simulation: ", step, "    From: ", Nsteps, "    Tmax: ", np.max(Temp), "   x: ", nv.xvec*1e+3 , " [mm]")
        #  ------------ Calculate Current Generated in Wire ------------ #
        V_current = TempPhysicalModels.CalculateCurrent(Flag_Current,Temp,1.0,dt)
        current1 = V_current[0]
        current2 = V_current[1]
        # -------------------------------------------------------------- #
      
        Time[step] = t
        Tmax[step] = np.max(Temp)

        Imax1[step] = current1
        Imax2[step] = current2

        FancyTemp = np.transpose(Temp) 

        # -------------- Update Position -------------- #
        if nv.WIRESCAN_Plane == "Horizontal":
            nv.WIRESCAN_wCposX += nv.WIRESCAN_wSpeed*dt
            V_Pos[step] = nv.xvec[0]
    
        elif nv.WIRESCAN_Plane == "Vertical":
            nv.WIRESCAN_wCposY += nv.WIRESCAN_wSpeed*dt
            V_Pos[step] = nv.yvec[0]
          
        TargetGeometry.CreateDetector(nv.DetType)

    nv.WireExp = TempPhysicalModels.LinearThermalExpansion(FancyTemp)

    return Time, Tmax, FancyTemp, Imax2, V_Pos