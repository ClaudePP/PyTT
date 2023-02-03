
'''
    This File Contains all the Variables Necessary to execute the codes. 
    This is done in this way so all the variables can be modified whenever necessary 
    by all the different functions. 

'''

import numpy as np
from Modules import ParticleBank as pb
from Modules import MaterialBank as mb

################################################
#            Physical Constants                #
################################################

ST = 5.6704E-8 		    # [J/sec m2 K4] Stefan-Bolthman constant for black-body radiation law
BZ = 1.38065e-23		# [J/K] Bolzmann constant
RH = 120.173			# [A/cm2 K2]	Richardson constant
h = 6.626e-34			# [m2kg/s] Planks constant
Pmass = 1.67E-27 		# [Kg] Mass of a proton
Emass = 9.1E-31 		# [Kg] Mass of an electron
ratio = Emass/Pmass	    # ratio
Na = 6.022E+23		    # Avogadro's Number
vluz = 2.998E+8		    # [m/s] Ligth velocity
pmass = 1.007316		# Proton relative mass
Qe = 1.60218e-19		# [C] elementary charge
re = 2.817940289e-15	# [m] Clasical electron radious
Na = 6.602e+23          # [part/mol] Avogadro's number
Um = 1.9944e-26         # [kg/uma] Mass of uma
Amu = 6.022e+26         # [uma/kg]

############################
ParticleInfo_dir = "asd.txt"

################################################
#        Define the beam  Parameters           #
################################################

RealInputFilename = "asd.txt"

BeamType = "BeamDescription/asd.txt"
ParticleProportionMatrix = np.array([])              # Matrix With Particle Porportion at each detector point. [posx,posy,prop]
Mat_BeamShape = []                                   # Matrix With beam shape (if not Gaussian)
Particle = pb.Particle('ParticleInfo/Default.txt')   # Name of particle the beam is made of. Loot in ParticleBank for available options. 
BEnergy = 0                                          # Beam Energy [MeV]
Intensity = 0                                        # Beam Intensity [A]
Nparticles = 0                                       # Number of particles per beam bunch.
sigx = 0                                             # sigma x of the beam [mm]
sigy = 0                                             # sigma y of the beam [mm]
tpulse = 0                                           # Length of the beam pulse [us]
frec = 0                                             # Frequency of the beam. Beam per second [Hz]
Npulses = 0                                          # Total number of pulses. 
x0 = 0.0                                             # Central position of the beam x [mm]
y0 = 0.0                                             # Central position of the beam y [mm]

################################################
#           Material Information               #
################################################

#Material = mb.Material('MaterialInfo/Default.txt')       # Name of Material detector is mada of. Look in Material Bank for available options.
Material = 0
enemat = 0                              # Energy deposited by Particle in Material. [MeV cm2/g]
Ele_enemat = 0                          # Energy Deposited by Particle's Electrons in Material [MeV cm2/g]

################################################
#         Geometry Information                 #
################################################

DetType = 'yxc'                     # Name of Detector Type. See TargetGeometry for available options. 

# Boolean variable for selecting detector type. 0: No 1: Yes
is_it_SEM = 0                        
is_it_WIRE = 0
is_it_FOIL = 0

# For the simulation, detectors will be simulated as a grid. 
xvec = np.asanyarray([0])           # Vector of x positions of each detor slice.
yvec = np.asanyarray([0])           # Vector of y postions of each detector slice. 

IntSurf = 0                         # [m2] Intensity surface: Surface of the detectors' slice particles see. 
eVol = 0                            # [m3] Detector's slice volume.
eSup = 0                            # [m2] Real surface of each detector slice.  

# -------   SEM  Grid Necessary Parameters  ----------

SEM_Plane = "asd"                   # Plane of beam to be measured: Horizontal, Vertical.
SEM_nWires = 0                      # Number of wires. 
SEM_wWidth = 0                      # Width of the wire. [m]
SEM_wLength = 0                     # Length of the wire [m]
SEM_wSep = 0                        # Separation between wires [m]
SEM_wRes = 0                        # Resolution of the wires: Each wire will be divided in slices. 
                                    # This is how big each slice will be. 
SEM_wSublim = np.asanyarray([0])    # Amount of sublimation suffered by the wire. [m]
                                    # it might be a vector. Each wire will srink depending 
                                    # on the maximum temperature reached by the wire

# -------   FOIL Necessary Parameters ------------------

FOIL_xwidth = 0                     # Length in x axis [m]
FOIL_nx = 0                         # Number divisions x axis
FOIL_ywidth = 0                     # Length in y axis [m]
FOIL_ny = 0                         # Number divisions y axis
FOIL_zwidth = 0                     # Length in z axis [m]


# -------    WIRE Necessary Parameters  ----------------------

WIRESCAN_Plane = "asd"               # Plane of beam to be measured: Horizontal, Vertical, Diagonal. 
WIRESCAN_wWidth = 0                  # Wire width [m]
WIRESCAN_wLength = 0.0               # Wire Length [m]
WIRESCAN_wRes = 0.0                  # Wire Resolution [m]
WIRESCAN_Type = 0.0                  # Type of Wire Scanner Moldel: 1 Continuous Beam. 2 Pulsated beam.
WIRESCAN_wSpeed = 0.0                # Wire velocity [m/s]
WIRESCAN_wCposX = 0.0                # Position of the center of the wire [m]
WIRESCAN_wCposY = 0.0                # Position of the center of the wire [m]
WIRESCAN_IniPos = 0.0                # Starting Scan Position [m]. Depending on Plane, it will be x or y. 
WIRESCAN_EndPos = 0.0                # Endins Scan Position [m]

################################################
#           Simulation Information             #
################################################

# ---------------- Temperature Simulation ---------------------- #
Flag_Temperature = 0
T0 = 0                              # Initial temperature of detector [K]
dtPulse = 0.0                       # Length Time Step Heating [s]
dtCooling = 0.0                     # Length Time step cooling [s]
EnableParameterVariation = 0        # Active: 1      Inactive: 0
RadiativeCooling = 0                # Active: 1      Inactive: 0
ThermoionicCooling = 0              # Active: 1      Inactive: 0
ConductiveCooling = 0               # Active: 1      Inactive: 0
SublimationCooling = 0              # Active: 1      Inactive: 0

# ---------------- Intensity Simulation ---------------------- #
Flag_Intensity = 0
Mu = 0.0
Eta = 0.0
BEp = 0.0
BEe = 0.0

################################################
#              Other Information               #
################################################

Mat_TempVapPress = []               # Matrix of Temperature vs vapor pressures. For sublimation calculation. 
WireExp = []                        # List of wire expansions. For SEM grid calculation of elongation. 

################################################
#              Output  Parameters              #
################################################

V_Time = []                         # Contains a list of all the time steps.
V_MaximumTemperature = []           # Contains a list of the maximum temperature reached by the detector at each instant of time. 
M_MaxTemp = []                      # Matrix, Maximum temperature of each wire along time. 
M_FancyTemperature = []             # Matrix, contains the maximum temperature reached by all the cells in the grid. 
V_Current1 = []                  # List of Currens along time for maximum current wire.  Without Thermoionic ems.
V_Current2 = []                  # List of Currens along time for maximum current wire.  With thermoionic ems.
M_Current = []                      # Matrix of currents along time. 
V_Emissivity = []                   # List of values of emissivity used. 
V_Pos = []                          # List of positions, used by wire scanners. 


CoolingImportance_Temp = []
CoolingImportance_Ems = []
CoolingImportance_Jth = []
CoolingImportance_Con = []
CoolingImportance_Sub = []

OutputFolderName = "Output/"


################################################
#                 Error Flags                  #
################################################

Flag1 = 0

Flag_QtInterface = 0


