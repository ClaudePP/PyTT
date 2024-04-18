
import numpy as np
from Modules import NecessaryVariables as nv
from Modules import TargetGeometry
import matplotlib.pyplot as plt

nv.WIRESCAN_Plane == "Horizontal"
nv.WIRESCAN_wLength = 0.05  # [m]
nv.WIRESCAN_wWidth = 30e-6  # [m]
nv.WIRESCAN_wRes = 30e-6 
nv.WIRESCAN_wCposX = 0.0
nv.WIRESCAN_wCposY = 0.0

TargetGeometry.CreateDetector('WIRESCAN')


rset=np.arange(0,0.0001,0.000001)
Uset=[]
for r in rset:
	Uset.append(TargetGeometry.BiasedWirePotential(5,r))
	print(r,Uset[-1])


plt.plot([r*1000 for r in rset],Uset,'o-')
plt.xlabel('distance [mm]')
plt.ylabel('potential [V]')
#plt.yscale('log')
#plt.ylim(0.1,10)
plt.show()