# ----------------------- Particle Banck ----------------------# 
#
# This file is very similar to the Material Banck module. It contains the information 
# of the particles conforming the beam. The Particle is a class, and we will create an 
# object of this class during the simulation. The properties of this particles are 
# properties of this calss. 
# 

import sys 

class Particle: 
    
    # This function will be called when creating an object of the class, it needs the 
    # input file with information of the particle properties. See the folder ParticleInfo for examples 
    # of particle file informations. 
    # 
    
    def __init__(self, ParticleFileName):
        
        h = open(ParticleFileName)

        d_PartInfo = {}
        cont = 0
        for l in h:
            if cont == 0:
                cont+=1
                continue
            else:
                d_PartInfo.update({l.split()[0] : l.split()[1]})

        self.name = d_PartInfo["Name:"]
        self.PartMass = float(d_PartInfo["ParticleMass:"])    # Particle mass  [Kg]
        self.Nprotons = float(d_PartInfo["Nprotons:"])        # Number of Protons.
        self.Nelectrons = float(d_PartInfo["Nelectrons:"])     # Number of electrons. 

        h.close()

      



    

    
