import numpy as np
import math
import sys
from Modules import MaterialBank as mb
from Modules import NecessaryVariables as nv

def CreateNiMatrix():
    def FindMatrixValue(x,y,Mat):
        if ((x < Mat[0][0]) or (x > Mat[-1][0])): val = 0
        elif ((y < Mat[0][1]) or (y > Mat[-1][1])): val = 0
        else:
            for i in range(0,len(Mat)):
                if ((x < Mat[i][0]) and (y < Mat[i][1])): break
            
            val = Mat[i][2]
        return val

    Nmat = np.zeros([len(nv.xvec),len(nv.yvec)])

    if nv.BeamType == "Gaussian":
        for k in range(len(nv.xvec)):
            for j in range(len(nv.yvec)):
                Ni = (1.0/(2*math.pi*nv.sigx*nv.sigy))*np.exp(-0.5*(((nv.xvec[k]-nv.x0)/nv.sigx)**2+((nv.yvec[j]-nv.y0)/nv.sigy)**2))
                Nmat[k,j] = Ni*1e-4
    else: 
        for i,x in enumerate(nv.xvec,0):
            for j,y in enumerate(nv.yvec,0):
                # first find corresponding element in BeamMatrix. 
                val = FindMatrixValue(x,y,nv.Mat_BeamShape)
                Nmat[i][j] = val*1e-4
    
    
    return Nmat


def NumberPartcles(Ntot):
    '''
    Gives you the number of particles reaching each part of the detector when the beam..
    '''
    Nmat = Ntot*nv.ParticleProportionMatrix
    
    return Nmat


def CalculateCurrent(Npart,Temperature,numberStepPulse,dt):
    '''
    Calculates the Current generated in each wire/foil.
    '''
    # PER INCIDENT PARTICLE 

    # Charge Deposition
    Qdep = nv.Particle.Nprotons*nv.Eta - nv.Particle.Nelectrons*nv.Mu

    # Secondary Emission
    Nv = nv.Na*nv.Material.rho/nv.Material.Am
    Ls = 1.0/(3.68e-17*Nv*nv.Material.Z**(1/3.))

    SEYp = 0.01*Ls*nv.enemat*1e+6*nv.Material.rho*(1.0+1.0/(1.0+(5.4*float(nv.BEnergy)/(nv.Particle.PartMass*nv.Amu))))
    SEYe = 0.01*Ls*nv.Ele_enemat*1e+6*nv.Material.rho

    Qse_p = nv.Particle.Nprotons*(1-nv.Eta)*SEYp + nv.Particle.Nprotons*nv.BEp*SEYp 
    Qse_e = nv.Particle.Nelectrons*(1-nv.Mu)*SEYe + nv.Particle.Nelectrons*nv.BEe*SEYe 
    
    # Total
    Q = Qdep + Qse_p + Qse_e

    # TOTAL CURRENT

    # Current Due to thermoionic emission. 
    thcurrent = nv.eSup*1e+4*nv.RH*Temperature**2*np.exp(-nv.Material.wfun*1.602e-19/(nv.BZ*Temperature))   # Current [A]

    # Number of particles reaching each spot in the geometry
    nparts = Npart*NumberPartcles(nv.Nparticles) / numberStepPulse
    
    # Calculate Current in each point of space due to Charge deposition and SEY. [A/m2]

    Super_Q = nparts*Q*nv.Qe/dt
    
    # Calculate Intensity.
        # Current 1: Without Thermoionic Emission. 
        # Current 2: With Thermoionic Emission. 

    if nv.DetType == "SEM":
        Surf = nv.SEM_wWidth*nv.SEM_wRes
        Current1 = []; Current2 = []
        if nv.SEM_Plane == "Vertical":
            Super_Q = Super_Q.copy().transpose()
        for k in range(0,len(Super_Q)):
            Current1 += [Surf*np.sum(Super_Q[k,:])]                           # Current  Without Thermoionic emission [A]
            Current2 += [Surf*np.sum(Super_Q[k,:])+np.sum(thcurrent[k,:])]      # Current With Thermoionic Emission [A] 
   
    elif nv.DetType == "FOIL":
        Surf = nv.FOIL_xwidth/nv.FOIL_nx * nv.FOIL_ywidth/nv.FOIL_ny
        Current1 = 0.0; Current2 = 0.0
        for k in range(0,len(Super_Q)):
            for j in range(0,len(Super_Q[k])):
                Current1 += Surf*Super_Q[k][j]
                Current2 += Surf*Super_Q[k][j]+Surf*thcurrent[k][j]

    
    elif nv.DetType == "WIRESCAN":
        Surf = nv.WIRESCAN_wRes*nv.WIRESCAN_wWidth
        Current1 = Surf*np.sum(Super_Q[0,:])                            # Current  Without Thermoionic emission [A]
        Current2 = Surf*np.sum(Super_Q[0,:])+np.sum(thcurrent[0,:])       # Current With Thermoionic Emission [A] 

    return Current1, Current2


def BeamHeating(Temperature, numberStepPulse):

    '''
    Calculates the beam heating of each space segment at a given instant of time.
     Arguments:
        Temperature: Temperature Matrix of the geometry .
        numberStepPulse: We assume constant temperature increase during beam pulse, so it
            is needed the number of steps the pulse is divided into.
            IMPORTANT! ---> If Wire Scanner Simulations  numberStepPulse = dt
    Output: Matrix of the increase of temperature in all the space. 
    '''

    if (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Type == 1):
        dt = numberStepPulse
        nparts = NumberPartcles(nv.Nparticles)*dt*nv.frec
    else:
        nparts = NumberPartcles(nv.Nparticles) / numberStepPulse

    dtemp =  nparts *  (nv.enemat+nv.Ele_enemat*nv.Particle.Nelectrons*nv.Mu)*1e+6*nv.Qe / nv.Material.CpT
    
    
    return np.asanyarray(dtemp)


def RadiativeCooling(dt, Temperature):

    '''
     Here radiative cooling is calculated
     Arguments:
        dt: [s] Is considered the duration of the cooling process. Normally it is defined by default but if
            if the conditions don't change too fast can be increased in order to make the simulation faster.
        Temperature: [K] Temperature Matrix.
    :return: temperature reduction (negative sign)
    '''

    cp = nv.Material.CpT
    eps = nv.Material.epsT
    dene = nv.eSup * nv.ST * eps * (Temperature ** 4 - (nv.T0 ** 4) * Temperature ** 0) * dt
    dtemp = -dene / (cp * nv.eVol * nv.Material.rho * 1e+6)

    return dtemp


def ThermoionicCooling(dt,Temperature):
    '''
    Here thermoionic cooling is calculated
    :param dt: [s] time duration of the cooling process.
    :param Temperature: [K] Temperature Matrix
    :return: temperature reduction (negative sign)
    '''

    thcurrent = nv.eSup* 1e+4 *nv.RH*Temperature**2*np.exp(-nv.Material.wfun*1.602e-19/(nv.BZ*Temperature))
    dene = (nv.Material.wfun*1.602e-19+(2*nv.BZ*Temperature))*thcurrent*dt/nv.Qe
    dtemp = -dene/(nv.Material.CpT*nv.eVol*nv.Material.rho*1e+6)

    return dtemp


def ConductiveCooling(dt, Temperature):
    '''
    Here Conductive cooling is calculated
    :param dt: [s] time duration of cooling process.
    :param Temperature: [K] Temperature matrix. 
    :return: Temperature variation. 
    '''
    dtemp = 0.0*Temperature

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
        elif (nv.DetType == "WIRESCAN") and (nv.WIRESCAN_Plane == "Vertical"):
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


def LinearThermalExpansion(Temperature):

    dL_vec = []
    for j in range(0,len(Temperature)):
        dL = 0.0
        for i in range(0,len(Temperature[j])):
            dL += nv.SEM_wRes*nv.Material.expcoeff*(Temperature[j][i]-nv.T0)

        dL_vec.append(dL)
    
    return dL_vec
