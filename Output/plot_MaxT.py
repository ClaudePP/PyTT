# plotting maximum temperature evolution
# 2023.11.15 mariousz.sapinski@psi.ch


import matplotlib.pyplot as plt
import pandas as pd

df = pd.read_csv("MaxTempVSTime.txt", skiprows=3, delim_whitespace=True, names=["em","time","maxT"])

plt.plot([x/1e6 for x in df["time"]],df["maxT"])
plt.xlabel("time [s]")
plt.ylabel("maximum temperature [K]")
plt.show()