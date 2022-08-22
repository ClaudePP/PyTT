'''
    This file contains the information about all the abailable particle types for 
    the beam. 
    They have been implemented in classes. 
'''

import sys 

class Particle: 

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

      



    

    
