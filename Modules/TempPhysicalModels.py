# ------------------ Physical Models ----------------------------------- # 
# 
# This module contains all the functions that have to do with the thermal physics of our simulations:
# Heating and cooling processes, but also
# calculating the amount of particles heating the wire and the current generation processes.
#


import numpy as np
import math
import sys
from Modules import MaterialBank as mb
from Modules import NecessaryVariables as nv
from scipy.integrate import dblquad  # 20240519, for beam profile integration
from scipy.constants import physical_constants


# global variable (constant)
elementary_charge = physical_constants["elementary charge"][0]
Stefan_Boltzmann = physical_constants["Stefan-Boltzmann constant"][0]




# ------------------------- Creating Particle Matrix ----------------------- # 

# The proportion of particles heating each point in space is stored in a matrix. 
# This matrix has thhe same size as the detector geometry, at the end we are only 
# interested in knowing how many particles heat each part of the detector. 
# This matrix is constructed either by reading the input file with the particle 
# distribution, or it is calculated from a gaussian function. 
#

def CreateNiMatrix():
    """
    this function is called multiple times
    every time it returns a beam distribution along the wire
    so it is a vector, not a matrix!    
    """
    def FindMatrixValue(x,y,Mat):
        if ((x < Mat[0][0]) or (x > Mat[-1][0])): val = 0
        elif ((y < Mat[0][1]) or (y > Mat[-1][1])): val = 0
        else:
            for i in range(0,len(Mat)):
                if ((x < Mat[i][0]) and (y < Mat[i][1])): break
            
            val = Mat[i][2]
        return val

    Nmat = np.zeros([len(nv.xvec),len(nv.yvec)])
    #print(len(nv.yvec))
    if nv.BeamType == "Gaussian":
        for k in range(len(nv.xvec)):
            for j in range(len(nv.yvec)):
                Ni = (1.0/(2*math.pi*nv.sigx*nv.sigy))*np.exp(-0.5*(((nv.xvec[k]-nv.x0)/nv.sigx)**2+((nv.yvec[j]-nv.y0)/nv.sigy)**2))
                # change ms: 20240405
                #Nmat[k,j] = Ni*1e-4  # sigx,y are in meters so Ni is in 1/m2 , why convert to 1/cm2?
                Nmat[k,j] = Ni  # sigx,y are in meters, Ni is 1/m2                
    else: 
        for i,x in enumerate(nv.xvec,0):
            for j,y in enumerate(nv.yvec,0):
                # first find corresponding element in BeamMatrix. 
                val = FindMatrixValue(x,y,nv.Mat_BeamShape)
                Nmat[i][j] = val
    
    # debugging
    if nv.Debug=="Beam":
        with open('Output/Nimatrix.txt','a') as fp:
            fp.write("Sum: "+str(Nmat.sum())+'\n')
            for iy in range(len(nv.yvec)):
                for ix in range(len(nv.xvec)):
                    fp.write(str(Nmat[ix][iy])+",")
                fp.write("\n")    
            #fp.write(str(Nmat[0][5]+','+str(Nmat.sum()))+'\n')
            #fp.write(str(Nmat.sum())+'\n')
    # this shows roughly corectness
    #exit()
    # to do: compute correction due to the binning!!!
    # to do: skew gaussian option 
    # unit of Nmat matrix is particles/m2
    return Nmat


# --------------------------------- Calculate Number of Particles ---------------------------- #

#
# Now we multiply the Matrix that was calculated in the previous funtion by the total number of particles. 
# This gives us the total number of particles heating each point in the detector. 
#

def NumberParticles(Ntot):
    '''
    Gives you the number of particles reaching each part of the detector 
    in fluence unit 1/m2.
    For a wire one dimension (wire movement direction) is always 1.
    This function is called twice per simulation step:
        1. BeamHeating
        2. CalculateCurrent
    '''
    Nmat = Ntot*nv.ParticleProportionMatrix
    
    return Nmat


# ms: new functions, calculates number of particles per second passing each segment of the wire 
# in case of gaussian beam
# define normalized 2D gaussian, with default beam parameters (center and sigmas)
def Gauss2D(x, y):
    '''
    2D normalized gaussian function

    Parameters
    ----------
    x : float
        horizontal position 
    y : float
        vertical position
    mx : float, optional
        beam center - horizontal coordinate. The default is nv.x0 = 0.
    my : float, optional
        beam center - vertical coordinate. The default is nv.y0 = 0.
    sx : float, optional
        horizontal sigma of the beam. The default is nv.sigx.
    sy : TYPE, optional
        vertical sigma of the beam. The default is nv.sigy.

    Returns
    -------
    TYPE
        Returns value of the normalized Gaussian distribution at given (x,y) coordinates.

    '''
    return 1. / (2. * np.pi * nv.sigx * nv.sigy) * np.exp(-((x - nv.x0)**2. / (2. * nv.sigx**2.) + (y - nv.y0)**2. / (2. * nv.sigy**2.)))
    

def ParticleRate():
    '''
    
    Returns
    -------
    float
        array of number of particles per second (npps) in each segment of the wire

    '''
    if nv.BeamType != "Gaussian":
        print('Cant use this function for nongaussian beam')
        exit(0)
    else:
        # nv.xvec and yvec are defined in Modules/TargetGeometry.py
        # npps is number of particles per second of Gaussian beam
        npps = np.zeros([len(nv.xvec),len(nv.yvec)])   # keep the same scheme as for old calculation
        if nv.DetType == "WIRESCAN":
            for k,xpos in enumerate(nv.xvec):
                #print("k=",k,' xpos = ',xpos)
                #exit(0)
                for j,ypos in enumerate(nv.yvec):
                    # wire moves horizontally and measures horizontal beam profile:
                        if nv.WIRESCAN_Plane =="Horizontal": 
                            x0=xpos-nv.WIRESCAN_wWidth/2
                            x1=xpos+nv.WIRESCAN_wWidth/2
                            y0=ypos-nv.WIRESCAN_wRes/2
                            y1=ypos+nv.WIRESCAN_wRes/2                    
                        elif nv.WIRESCAN_Plane =="Vertical": 
                            x0=xpos-nv.WIRESCAN_wRes/2
                            x1=xpos+nv.WIRESCAN_wRes/2
                            y0=ypos-nv.WIRESCAN_wWidth/2
                            y1=ypos+nv.WIRESCAN_wWidth/2   
                        #print('x0,x1 = ',x0,x1)    # Debug  
                        #print('y0,y1 = ',y0,y1)    # Debug
                        # x and y are exchanged, this gives proper orientation of gaussian.
                        npps[k,j],err=dblquad(Gauss2D,y0,y1,lambda x: x0,lambda x: x1,epsrel = 1e-10, epsabs = 0)
                        #print("nv.xvec=",nv.xvec)  # debug
                        #print("nv.yvec=",nv.yvec)  # debug                        
                        #print('npps [k,j] = ',npps[k,j]) # debugging
        if nv.DetType == "SEM":
            # add number of wires!
            for k,xpos in enumerate(nv.xvec):
                #print("Debug TempPhysicalModel k=",k,' xpos = ',xpos)
                #exit(0)
                for j,ypos in enumerate(nv.yvec):
                    # wire moves horizontally and measures horizontal beam profile:
                        if nv.SEM_Plane =="Horizontal": 
                            x0=xpos-nv.SEM_wWidth/2
                            x1=xpos+nv.SEM_wWidth/2
                            y0=ypos-nv.SEM_wRes/2
                            y1=ypos+nv.SEM_wRes/2                    
                        elif nv.SEM_Plane =="Vertical": 
                            x0=xpos-nv.SEM_wRes/2
                            x1=xpos+nv.SEM_wRes/2
                            y0=ypos-nv.SEM_wWidth/2
                            y1=ypos+nv.SEM_wWidth/2   
                        #print('x0,x1 = ',x0,x1)    # Debug
                        #print('y0,y1 = ',y0,y1)    # Debug
                        # x and y are exchanged, this gives proper orientation of gaussian.
                        npps[k,j],err=dblquad(Gauss2D,y0,y1,lambda x: x0,lambda x: x1,epsrel = 1e-10, epsabs = 0)
                        #print("nv.xvec=",nv.xvec)   # Debug
                        #print("nv.yvec=",nv.yvec)   # Debug                     
        if nv.Debug=="Beam":
            print("TempPhysicalModels:ParticleRate npps=",npps)    # Debug, controlled January 2025, OK, it makes sense!
            # npps is an array with dimentsion [nwires][nbins]
        # nv.Intensity is in Amps    
    return npps*nv.Intensity/elementary_charge # assuming particles with charge 1 (protons)




#------------------------- Current Calculation -------------------------------------------- # 

# 
# This function calculates the secondary electron and thermionic currents 
# generated in the detector at each time step. 
# Npart is just 1.0 or 0.0

def CalculateCurrent(Npart,Temperature,numberStepPulse,dt):
   
    # First we calculate how much charge each incident generates in the material. 
    # Qtot = Qdep + Qse + Qje
    # Total charge accounts for the deposited charge, the secondary emission charge and the thermionic charge. 
    # 
    
    # 1) Deposited charge, depends on the amount of electrons and protons that stay in the detector material. 
    #
    
    Qdep = nv.Particle.Nprotons*nv.Eta - nv.Particle.Nelectrons*nv.Mu

    # 2) Secondary Emission, has to do with how many electrons leave the material per incident particle. 
    # It is characterized by the SEY, secondary emission yield. We have the SEY generated by the incident 
    # protons and the SEY generated by the incident electrons. 
    # 
    
    Nv = nv.Na*nv.Material.rho/nv.Material.Am   # [atoms/cm3] atomic density
    # Ls is an effective penetration distance of secondary electrons 
    # 3.68e-17 is a constant with unit of [cm2], see for instance D.Kramer thesis page 29-30
    Ls = 1.0/(3.68e-17*Nv*nv.Material.Z**(1/3.))  # [cm]

    # ms, 20230808:    
    #SEYp = 0.01*Ls*nv.enemat*1e+6*nv.Material.rho*(1.0+1.0/(1.0+(5.4*float(nv.BEnergy)/(nv.Particle.PartMass*nv.Amu))))
    # nv.Material.rho is in [g/cm3]
    # nv.enemat is dE/dx in [MeV*cm2/g] and it needs to be converted to eV/cm, hence factor 1e+6
    if nv.Debug=="Edep":
        print("Debug   TempPhysicalModels:CalculateCurrent:nv.enemat: ",nv.enemat)
    dEdxel=nv.enemat*1e+6*nv.Material.rho   # [eV/cm] electronic energy loss
    BEnergy_eV = float(nv.BEnergy)*1e6
    pmassamu = 1.00727647
    # units are:
        # 0.01 - [1/eV]
        # Ls - [cm]
        # dEdxel - [eV/cm]
        # 5.4e-6 - [amu/eV]
    SEYp = 0.01*Ls*dEdxel*(1.0+1.0/(1.0+(5.4e-6*BEnergy_eV/pmassamu)))
    if nv.Debug=='SEY':
        print("Debug   TempPhysicalModels:SEY protons = ",SEYp)
        #exit(0)
    nv.S_SEYp=SEYp # to be saved in output
    # BTW delta electrons lead to Secondary electron emission!
    SEYe = 0.01*Ls*nv.Ele_enemat*1e+6*nv.Material.rho

    # ms, 20230808:
    #Qse_p = nv.Particle.Nprotons*(1-nv.Eta)*SEYp + nv.Particle.Nprotons*nv.BEp*SEYp 
    # nv.Eta - fraction of protons stopped in the material (thin targets: 0)
    # nv.BSp - fraction of backscatterted protons (they also cross the wire surface twice)
    # Idea: plus additional factor when they "scratch" edge of the wire?
    # BTW this equation takes into account the fact, that secondary electrons are produced on two surfaces: entry and exit.
    # remark by Araceli: the entry surface is taken into account by the last term, 
    # the exit surface is considered in the first term; when eta=0 protons are not stopped and secondary electrons is emitted
    Qse_p = nv.Particle.Nprotons*(1-nv.Eta)*SEYp + nv.Particle.Nprotons*nv.BEp*SEYp + nv.Particle.Nprotons*SEYp   
    Qse_e = nv.Particle.Nelectrons*(1-nv.Mu)*SEYe + nv.Particle.Nelectrons*nv.BEe*SEYe 
    
    
    # Here we add upp the aforecalculated charges, deposited charge and SE charge. 
    # Q is number of charges generated per impacting beam particle
    
    Q = Qdep + Qse_p + Qse_e  # 2024.04.05 - check if here should be minus Qse_e
    if nv.Debug=='SEY':
        print("Debug   TempPhysicalModels:Q = ",Q)
        #exit(0)


    # Current Due to thermionic emission. This term highly depends on the temperature of the detector at each point. 
    #
 
    # not needed, remove:
    #if nv.EnableParameterVariation:
    #    nv.Material.wf = nv.Material.GetWf(Temperature)            
    #    print(Temperature)
    #    print(nv.Material.wf)
    #else: 
    #    nv.Material.wf = nv.Material.GetWf(300*Temperature**0)
 
    
    # eSup is vector of surfaces of each bin of the wire [m2]
    # RH is Richardson constant 120.173	 [A/cm2 K2] - factor 1e4 converts it to to [A/m2 K2]   
    # wfun is material work function [eV]
    if nv.Debug=="Thermionic":
        print("Debug   TempPhysicalModels:CalculateCurrent:thcurrent:nv.Material.wfun: ",nv.Material.wfun)
        print("Debug   TempPhysicalModels:CalculateCurrent:thcurrent:nv.RH: ",nv.RH)
        print("Debug   TempPhysicalModels:CalculateCurrent:thcurrent:nv.eSup: ",nv.eSup)
        print("Debug   TempPhysicalModels:CalculateCurrent:thcurrent:nv.Temperature: ",Temperature)
    #thcurrent = nv.eSup*1e+4*nv.RH*Temperature**2*np.exp(-nv.Material.wfun*nv.Qe/(nv.BZ*Temperature))   # Current [A]
    thcurrent = nv.eSup*1e4*nv.RH*Temperature**2*np.exp(-nv.Material.wfun*nv.Qe/(nv.BZ*Temperature))   # Current [A]

    
    #
    # Number of particles reaching each spot in the geometry, it is used to calculate the total current. 
    
    # nparts = Npart*NumberParticles(nv.Nparticles) / numberStepPulse
    # new, 2024.05.20: 
    nparts = ParticleRate()*dt   # this only works for WireScanner for now
        
    # Calculate Current in each point of space due to Charge deposition and SEY. [A/m2]
    # Here we transform the charge calculated before to current. 

    Super_Q = nparts*Q*nv.Qe/dt   # [A]
    if (nv.Debug=="Electric") or (nv.Debug=="SEY"):
        print("TempPhysicalModels:CalculateCurrent Super_Q: ",Super_Q)
        print("TempPhysicalModels:CalculateCurrent len(Super_Q): ",len(Super_Q))
        for k in range(0,len(Super_Q)):
            print(np.sum(Super_Q[k,:]))
    
    
    # At this stage we have a matrix with the amount of current generated in each point in space, 
    # but we are interested in how much total current is generated in the detector, so we need to add up 
    # all this currents. 
    # This is slightly different for every type of detector, for example in a wire scanner we will integrate 
    # along the single wire lenght. For a SEM grid we have to integrate each wire sepparatelly in case there are 
    # more than one. And in the foil we need to integrate both spatial coordinates. 
    # We have sepparated the current in two different ones. 
    # Current1 does not account for thermionic emission. Current2 accounts for thermionic emission. This might 
    # give us information about the thermal contribution. 
    #
    

    if nv.DetType == "SEM":
        Surf = nv.SEM_wWidth*nv.SEM_wRes
        Current1 = []; Current2 = []
        if nv.SEM_Plane == "Vertical":
            Super_Q = Super_Q.copy().transpose()
        for k in range(0,len(Super_Q)):  # this is a loop over wires in fact
            #Current1 += [Surf*np.sum(Super_Q[k,:])]                           # Current  Without Thermionic emission [A]
            #Current2 += [Surf*np.sum(Super_Q[k,:])+np.sum(thcurrent[k,:])]      # Current With Thermionic Emission [A] 
            # 20250131: , remove Surf?  ****
            #Current1 += [Surf*np.sum(Super_Q[k,:])]                           # SEM current [A]
            #Current2 += [Surf*np.sum(thcurrent[k,:])]                         # Thermionic current [A] 
            Current1 += [np.sum(Super_Q[k,:])]                           # SEM current [A]
            Current2 += [np.sum(thcurrent[k,:])]                         # Thermionic current [A] 

        print("Current1: ", Current1)  # debug array of np.floats with dimension = number of wires
        print("Current2: ", Current2)  # debug array of np.floats with dimension = number of wires

    
    elif nv.DetType == "FOIL":
        Surf = nv.FOIL_xwidth/nv.FOIL_nx * nv.FOIL_ywidth/nv.FOIL_ny
        Current1 = 0.0; Current2 = 0.0
        for k in range(0,len(Super_Q)):
            for j in range(0,len(Super_Q[k])):
                Current1 += Surf*Super_Q[k][j]  # remove Surf
                Current2 += Surf*Super_Q[k][j]+Surf*thcurrent[k][j]

    
    elif nv.DetType == "WIRESCAN":
        # ms: 20230808
        # following Manon's suggestion:
        # Surf = nv.WIRESCAN_wRes*nv.WIRESCAN_wWidth*1e4            # 1e4? [cm2]?
        Surf = nv.WIRESCAN_wRes*nv.WIRESCAN_wWidth                  # [m2] - number of particles crossing the wire in each bin
        Current1 = np.sum(Super_Q[0,:])                             # [A] Current  Without Thermionic emission
        Current2 = np.sum(thcurrent[0,:])                           # [A] Thermionic Current
        #Current1 = np.sum(Super_Q[0,:])                            # Current  Without Thermionic emission [A]
        #Current2 = np.sum(thcurrent[0,:])                          # Thermionic Emission current [A] 
        # add Current3 delta electrons

    # ms: 20240418, initial coding, not validated
    elif nv.DetType == "SPLITTER":    # splitter can also be aperture foil
        Surf = nv.SPLITTER_wRes*nv.SPLITTER_wWidth                  # [m2] - number of particles crossing the wire in each bin
        Current1 = Surf*np.sum(Super_Q[0,:])                        # Current  Without Thermionic emission [A]
        Current2 = np.sum(thcurrent[0,:]) 

    return Current1, Current2

# -------------------------------------------- BEAM HEATING -------------------------------- # 
#
# Calculates the beam heating of each space segment at a given instant of time. We are considering that 
# the heating will be always the same during the beam pulse.
# 2024.03.29 - add sum of the number of particles to the output

def BeamHeating(Temperature, numberStepPulse):
    
    # Determine the number of particles nparts
    if (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Type == 1):
        dt = numberStepPulse
        # ms: 20230808, after Manon:
        #nparts = NumberParticles(nv.Nparticles)*dt*nv.frec   # synchrotron beam
        # this below should be from the current
        #nparts = NumberParticles(nv.Nparticles)               # cyclotron beam
        # ms: 20240519, any beam:
        nparts = ParticleRate()*dt     
    else:
        #nparts = NumberParticles(nv.Nparticles) / numberStepPulse # 2correct
        dt = nv.dtPulse        
        nparts = ParticleRate()*dt     
        #print("TempPhysicalModels:BeamHeating nparts = ", nparts, nparts.sum()," dt = ",dt) # debug (Beam)
        
    if nv.Debug=="Beam":       
        print(f"[DEBUG] TempPhysicalModels BeamHeating: nparts = {nparts}, sum = {nparts.sum()}, dt = {dt}")
    

    # Debugging to output file
    if nv.Debug=='Beam':
        fout='Output/beam_profile.txt'
        with open(fout,'a') as fo:
            #fo.write("beam matrix dimensions: "+str(len(nv.xvec))+','+str(len(nv.yvec))+'\n')
            #fo.write("simstep: "+str(len(nv.M_MaxTemp))+'\n')
            fo.write("simstep: "+str(nv.stepcount)+'\n')            
            for k in range(len(nv.xvec)):
                for j in range(len(nv.yvec)):
                    fo.write(str(k)+','+str(j)+','+str(nparts[k][j])+'\n')


    # Energy deposition and temperature increase calculation:
    # 2024.04.11: added 1e-4 because nparts in in 1/m2 now
    # nv.Material.CpT [J/(gK)]
    #dtemp =  nparts *  1e-4 * (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e+6*nv.Qe / nv.Material.CpT
    # 2024.05.19: removed 1e-4 because nparts is just number of particles
    # 1e6 converts MeV to eV, and nv.Qe converts them to Joules
    # nv.enemat is in MeV*g/cm2, CpT is in [K/g*J] 
    energy_factor = (nv.enemat + nv.Ele_enemat * nv.Particle.Nelectrons * nv.Mu) * 1e6 * elementary_charge
    material_cp = nv.Material.CpT  # Specific heat capacity [J/(gK)]
    if nv.DetType == "WIRESCAN":
        #nptsum=0
        #dtemp =  nparts * (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e+6*elementary_charge / (nv.WIRESCAN_wWidth*nv.WIRESCAN_wRes*nv.Material.CpT)
        dene = nparts * energy_factor * nv.Material.rho * nv.WIRESCAN_wWidth * 100  # J
        dtemp = nparts * energy_factor / (nv.WIRESCAN_wWidth * nv.WIRESCAN_wRes * material_cp)
        #dene = nparts * (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e6*elementary_charge*nv.Material.rho*nv.WIRESCAN_wWidth*100  # J
        nptsum=nparts.sum()
    if nv.DetType == "SEM":
        dtemp = nparts * energy_factor / (nv.SEM_wWidth * nv.SEM_wRes * material_cp)
        dene = nparts * energy_factor * nv.Material.rho * nv.SEM_wWidth * 100  # J
        nptsum = [np.sum(nparts[k, :]) for k in range(nv.SEM_nWires)]

        # rm:        
        #nptsum=[]
        #dtemp =  nparts * (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e+6*nv.Qe / (nv.SEM_wWidth*nv.SEM_wRes*nv.Material.CpT)
        #dene = nparts * (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e6*nv.Qe*nv.Material.rho*nv.SEM_wWidth*100  # J
        #for k in range(0,nv.SEM_nWires):
        #    nptsum.append(np.sum(nparts[k,:]))
        
    # Debugging output
    if nv.Debug=="Edep":
        print(f"[DEBUG] TempPhysicalModels:BeamHeating: dene = {dene}")

    # 1e6 is to convert g/cm3 to g/m3 
    dtemp = dene/(material_cp * nv.eVol * nv.Material.rho * 1e6)

    #print("Debug TempPhysicalModels:BeamHeating:nv.Material.CpT ",nv.Material.CpT)
    # As an output we obtain the temperature variation for each point.    
    # 2024.03.30: include also total number of particles in output, do we need it? We have nv.Nparticles
    # 2024.05.20: include energy in the output
    return np.asanyarray(dtemp),np.asarray(dene),nptsum

# ---------------------------------- Radiative Cooling -------------------------- # 

def RadiativeCooling(dt, Temperature):

    '''
     Here radiative cooling is calculated
     Arguments:
        dt: [s] Is considered the duration of the cooling process. Normally it is defined by default but if
            if the conditions don't change too fast can be increased in order to make the simulation faster.
        Temperature: [K] Temperature Matrix.
    :return: temperature reduction (negative sign)
    '''

    cp = nv.Material.CpT     # [J/(gK)]
    eps = nv.Material.epsT   # emmissivity []
    #Stefan-Boltzmann constant
    # nv.ST - Stefan-Boltzmann [J/sec m2 K4]
    dene = nv.eSup * Stefan_Boltzmann * eps * (Temperature ** 4 - (nv.T0 ** 4) * Temperature ** 0) * dt
    dtemp = -dene / (cp * nv.eVol * nv.Material.rho * 1e+6)

    return dtemp,dene

# ------------------------------- Thermionic Cooling --------------------------------------- # 
# 
def ThermionicCooling(dt,Temperature):
    '''
    Here thermionic cooling is calculated
    :param dt: [s] time duration of the cooling process.
    :param Temperature: [K] Temperature Matrix
    :return: temperature reduction (negative sign)
    '''
    
    # thcurrent should be a single function! It is defined twice.
    # 1.602e-19 - use Qe 
    #thcurrent = nv.eSup*1e4*nv.RH*Temperature**2*np.exp(-nv.Material.wfun*1.602e-19/(nv.BZ*Temperature))
    thcurrent = nv.eSup*1e4*nv.RH*Temperature**2*np.exp(-nv.Material.wfun*nv.Qe/(nv.BZ*Temperature))   # thermionic current [A]
    # 2024.05.20:
    #dene = (nv.Material.wfun*nv.Qe+(2*nv.BZ*Temperature))*thcurrent*dt/nv.Qe
    dene = (nv.Material.wfun*nv.Qe+(2*nv.BZ*Temperature))*thcurrent*dt/nv.Qe
    dtemp = -dene/(nv.Material.CpT*nv.eVol*nv.Material.rho*1e+6)

    return dtemp, dene


# ------------------------------- Thermionic Electron Spectra ----------------------------- # 
# 
def ThermionicElectronSpectra(Temp,Eks):
    '''

    Parameters
    ----------
    Temp : float
        target temperature [K]
    Eks : array of floats    
        array of energies [eV]
    Returns
    -------
    fracs: array of floats
        relative intensities for each energy in Eks
    '''
    kBeV=nv.BZ/nv.Qe
    fracs = (Eks/pow(2*kBeV*Temp,2))*np.exp(-Eks/(2*kBeV*Temp))
    
    return fracs


# ------------------------------- Electron velocity from kinetic energy -------------------- # 
# 
def ElectronVelocity(Eks):
    '''
    Computes electron velocity [m/s] from kinetic energy [eV].
    !!! Non-relativistic !!!

    Parameters
    ----------
    Eks : float
        electron kinetic energy [eV]

    Returns
    -------
    ve : float
        electron velocity [m/s]

    '''
    EJoule=Eks*nv.Qe  # [J] convert eV to Joules
    ve = np.sqrt(2*EJoule/nv.Emass)
    return ve



# --------------------------------- Conductive cooling --------------------------------- # 

def ConductiveCooling(dt, Temperature):
    '''
    Here Conductive cooling is calculated
    :param dt: [s] time duration of cooling process.
    :param Temperature: [K] Temperature matrix. 
    :return: Temperature variation. 
    '''
    
    dtemp = 0.0*Temperature

    # --------------- This is based on the FTCS method --------------- $
    
    def Calculate_dTemp(posvec1,posvec2,Temperature,dt,dx):
        Tbott = 300.0; Tupp = 300.0 
        for k in range(0,len(posvec1)):
            for j in range(0, len(posvec2)):
                Tj = Temperature[k][j]
                alpha = nv.Material.con[k][j]/(nv.Material.rho*nv.Material.CpT[k][j]*1e+6)
                r = alpha*dt/dx**2
                if j == 0:
                    Tjp1 = Temperature[k][j+1]
                    Tjm1 = Tbott
                elif j == len(posvec2)-1:
                    Tjp1 = Tupp
                    Tjm1 = Temperature[k][j-1]
                else:
                    Tjp1 = Temperature[k][j+1]
                    Tjm1 = Temperature[k][j-1]

                dtemp[k][j] = r*(Tjp1-2*Tj+Tjm1)
        return dtemp

    # Two possible cases, 
    #       for SEM Grids and Wire scanners only conduction in one direction. 
    #       for FOILs conduction in direction x and y. 

    if nv.DetType != "FOIL":
        #-----------------------------------------------#
        # posvec1: Position of Wire Center.             #
        # posvec2: Central position of wire fragment.   #
        #-----------------------------------------------#
        if (nv.DetType == "SEM") and (nv.SEM_Plane == "Horizontal"):
            posvec1 = nv.xvec.copy()
            posvec2 = nv.yvec.copy()
            dx = nv.SEM_wRes
            Calculate_dTemp(posvec1,posvec2,Temperature,dt,dx)
        # ms: 20230808, following Manon
        #elif (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Plane == "Vertical"):
        elif (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Plane == "Horizontal"):    
            posvec1 = nv.xvec.copy()
            posvec2 = nv.yvec.copy()
            dx = nv.WIRESCAN_wRes
            Calculate_dTemp(posvec1,posvec2,Temperature,dt,dx)
        elif (nv.DetType == "SEM") and (nv.SEM_Plane == "Horizontal"):
            posvec1 = nv.yvec.copy()
            posvec2 = nv.xvec.copy()
            dx = nv.SEM_wRes
            Temperature = Temperature.transpose()
            nv.Material.con = nv.Material.con.transpose()
            nv.Material.CpT = nv.Material.CpT.transpose()
            dtemp = dtemp.transpose()
            # ---------------------------------------------- #
            dtemp = Calculate_dTemp(posvec1,posvec2,Temperature,dt,dx)
            # ---------------------------------------------- #
            Temperature = Temperature.transpose()
            nv.Material.con = nv.Material.con.transpose()
            nv.Material.CpT = nv.Material.CpT.transpose()
            dtemp = dtemp.transpose()
        elif (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Plane == "Vertical"):
            posvec1 = nv.yvec.copy()
            posvec2 = nv.xvec.copy()
            dx = nv.SEM_wRes
            Temperature = Temperature.transpose()
            nv.Material.con = nv.Material.con.transpose()
            nv.Material.CpT = nv.Material.CpT.transpose()
            dtemp = dtemp.transpose()
            # ---------------------------------------------- #
            dtemp =  Calculate_dTemp(posvec1,posvec2,Temperature,dt,dx)
            # ---------------------------------------------- #
            Temperature = Temperature.transpose()
            nv.Material.con = nv.Material.con.transpose()
            nv.Material.CpT = nv.Material.CpT.transpose()
            dtemp = dtemp.transpose()

        
    else:
        # - Calculate conductive cooling for FOILS. By solving 2d Heat equation - #  
        #       Tij = T[i][j]; Tim1 = T[i-1][j]; Tip1 = T[i+1][j]                 #
        #                      Tjm1 = T[i][j-1]; Tjp1 = T[i][j+1]                 #
        # ----------------------------------------------------------------------- #

        dx = nv.FOIL_xwidth/nv.FOIL_nx; dy = nv.FOIL_ywidth/nv.FOIL_ny

        Tbott = 300.0; Tupp = 300.0 ; Tright = 300.0; Tleft = 300.0
        for i in range(0,len(nv.xvec)):
            for j in range(0,len(nv.yvec)):
                Tij = Temperature[i][j]
                alpha = nv.Material.con[i][j]/(nv.Material.rho*nv.Material.CpT[i][j] * 1e+6)
                rx = alpha*dt/dx**2; ry = alpha*dt/dy**2
                
                if (i == 0) and (j == 0): 
                    Tim1 = Tleft; Tjm1 = Tbott
                    Tip1 = Temperature[i+1][j]; Tjp1 = Temperature[i][j+1]
                elif (i == len(nv.xvec)-1) and (j == len(nv.yvec)-1): 
                    Tim1 = Temperature[i-1][j]; Tjm1 = Temperature[i][j-1]
                    Tip1 = Tright; Tjp1 = Tupp

                else:
                    try: Tim1 = Temperature[i-1][j]; Tip1 = Temperature[i+1][j] 
                    except IndexError:
                        if (i == 0): Tim1 = Tleft
                        elif (i == len(nv.xvec)-1): Tip1 = Tright
                    try: Tjm1 = Temperature[i][j-1]; Tjp1 = Temperature[i][j+1]
                    except IndexError:
                        if (j == 0): Tjm1 = Tbott
                        elif (j == len(nv.yvec)-1): Tjp1 = Tupp
                dtemp[i][j] = ry*(Tjp1-2*Tij+Tjm1)+rx*(Tip1-2*Tij+Tim1)
    return dtemp

# ------------------------------- Sublimation Cooling ---------------------------------- # 

def SublimationCooling(dt, Temperature):
    '''
        This function calculates the sublimation rate of the material due to temperature. 
        The model is based on " Scientific fundatios of Vacuum techniques", S. Dushman.  Eq. 10.09.
        :return: Temperature reduction (negative sign)
    '''
    dtemp = Temperature * 0.0
    for i in range(0,len(Temperature)):
        for j in range(0,len(Temperature[i])):
            Tij = Temperature[i][j]
            Log_W = nv.Material.Sublimation_C1 - 1/2.*np.log(Tij) - nv.Material.Sublimation_C2/Tij
            sub_W = math.pow(10.0, Log_W)                           # [g / cm2 / s]   material sublimation rate
            sub_Wdt = sub_W * dt * 10000                            # [g / cm2 / s] -> [g / m2]
            Wdt = sub_Wdt*nv.eSup * nv.Material.Am                  # [J]     
            dtemp[i][j] = -Wdt*nv.Material.HT[i][j] / (nv.Material.CpT[i][j] * nv.eVol * 1e+6 * nv.Material.rho )            # [ K ]

    return dtemp

# ------------------------  Linear thermal expanison --------------------------------- 

# This attempts to calculate the thermal expanison of the detecotrs, but it is not currently being used.  

def LinearThermalExpansion(Temperature):

    dL_vec = []
    for j in range(0,len(Temperature)):
        dL = 0.0
        for i in range(0,len(Temperature[j])):
            dL += nv.SEM_wRes*nv.Material.expcoeff*(Temperature[j][i]-nv.T0)

        dL_vec.append(dL)
    
    return dL_vec
