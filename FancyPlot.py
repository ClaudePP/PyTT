import numpy as np
import matplotlib.pyplot as plt

def readfile(filename):
    f = open(filename)
    Time, Temp = [],[]
    for k,l in enumerate(f,start=0): 
        if k > 2:
            Time += [float(l.split(" ")[3])*1e-6]
            Temp += [float(l.split(" ")[-1].split("\n")[0])]
            

    return Time, Temp
# ----------------------------------------- #

filename = "Output/MaxTempVSTime.txt"
Time, Temp = readfile(filename)

plt.figure(2,figsize=(8,6))
plt.grid()
plt.plot(Time,Temp,color="darkred",lw=2)
plt.title("Max. Temperature: CHROMOX",fontsize=14)
plt.xlabel("Time [s]",fontsize=14)
plt.ylabel("Temperature [K]",fontsize=14)
plt.ylim([300,600])
plt.xticks(fontsize=14); plt.yticks(fontsize=14)

plt.show()