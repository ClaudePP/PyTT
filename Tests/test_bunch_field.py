
import numpy as np
from Modules import NecessaryVariables as nv
from Modules import TargetGeometry
from Modules import BunchField
import matplotlib.pyplot as plt

nv.WIRESCAN_Plane == "Horizontal"
nv.WIRESCAN_wLength = 0.05  # [m]
nv.WIRESCAN_wWidth = 30e-6  # [m]
nv.WIRESCAN_wRes = 30e-6 
nv.WIRESCAN_wCposX = 0.0
nv.WIRESCAN_wCposY = 0.0

TargetGeometry.CreateDetector('WIRESCAN')

# wire, bias 30V
#rset=np.arange(-0.1,0.1,0.00001)
rset=np.arange(0.0,0.1,0.00001)

Uset=[]

# bunch
nv.sigx = 1.6e-3 # [m]
sigmaz=1e-2         # [m]
bcharge=1.0e8   # [nparts] 2.5e8
posz=0.0
Ebunch=[]

# calculate wire potential and bunch field
for r in rset:
	Uset.append(-1*TargetGeometry.BiasedWirePotential(60,abs(r+0.001)))
	Ebunch.append(BunchField.SymmetricBunchField(abs(r), posz, sigmaz, bcharge))
	print(r,Uset[-1],Ebunch[-1])

# calculate bunch potential
dr=rset[1]-rset[0]
print("dr [m] = ",dr)
Ubunch=[]
for idx,ef in enumerate(Ebunch):
	Ubunch.append(-1*dr*np.cumsum(Ebunch[idx:])[-1])
#exit()


# negative r:
rsetn=np.arange(-0.1,0.0,0.00001)
Usetn=[]
Ebunchn=[]
for r in rsetn:
	Usetn.append(-1*TargetGeometry.BiasedWirePotential(60,abs(r+0.001)))
	Ebunchn.append(BunchField.SymmetricBunchField(abs(r), posz, sigmaz, bcharge))

Ubunchn=[]
for idx,ef in enumerate(Ebunchn):
	if idx>0:
		Ubunchn.append(-1*dr*np.cumsum(Ebunchn[:idx])[-1])

Ubunchn.append(Ubunchn[-1])


# sum
rsett=[]
Utot=[]
for idx,r in enumerate(rsetn):
	rsett.append(r)
	Utot.append(Ubunchn[idx]+Usetn[idx])
for idx,r in enumerate(rset):
	rsett.append(r)
	Utot.append(Ubunch[idx]+Uset[idx])


plt.plot([(r*1000) for r in rset],Uset,'b-',label='wire potential') # wire potential
plt.plot([(r*1000) for r in rsetn],Usetn,'b-') # wire potential
#plt.plot([r*1000 for r in rset],Ebunch,'o-') # bunch field
plt.plot([r*1000 for r in rset],Ubunch,'go-',label='bunch potential') # bunch potential
plt.plot([r*1000 for r in rsetn],Ubunchn,'go-') # bunch potential
#plt.plot([r*1000 for r in rsett],Utot,'r',label="sum")
plt.xlabel('distance [mm]')
plt.ylabel('potential [V]')
plt.legend()
plt.xlim(-9,5)
plt.ylim(-65,-25)
plt.grid()
plt.show()