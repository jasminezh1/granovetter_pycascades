import numpy as np
import matplotlib
matplotlib.use('Agg') #otherwise use 'pdf' instead of 'Agg'
import matplotlib.pyplot as plt
from matplotlib import colors
import re
import glob
import os
import math

import seaborn as sns
sns.set(font_scale=1.5)
sns.set_style("whitegrid")
sns.despine()
import scipy.io


def co2_gmt(co2):
    alpha = 4.841 # see IPCC, AR6, chapter 6, table A 6.2
    beta = 0.0906
    co2_preind = 280
    #factor = 0.266 #where does this come from?
    gmt = alpha*math.log(co2/co2_preind) + beta*(np.sqrt(co2)-np.sqrt(co2_preind))
    return gmt


#read in matlab data file
mat = scipy.io.loadmat("Scenarios.mat")


#CO2 values
x_MF = mat['x_MF'].flatten()
x_WF = mat['x_WF'].flatten()
x_SF = mat['x_SF'].flatten()




#Conversion to GMT
GMT_MF = []
GMT_WF = []
GMT_SF = []
for i in range(0, len(x_MF)):
    gmt_val_MF = co2_gmt(x_MF[i])
    gmt_val_WF = co2_gmt(x_WF[i])
    gmt_val_SF = co2_gmt(x_SF[i])
    GMT_MF.append(gmt_val_MF)
    GMT_WF.append(gmt_val_WF)
    GMT_SF.append(gmt_val_SF)
print(GMT_MF[0])
GMT_MF = np.array(GMT_MF)/GMT_MF[0] #the 1/GMT_MF[0] is the rescaling factor to set 408ppm to 1.0°C above pre-industrial, 0.266 does not make sense to get reasonable GMT-values
GMT_WF = np.array(GMT_WF)/GMT_WF[0]
GMT_SF = np.array(GMT_SF)/GMT_SF[0]


np.savetxt("GMT/gmt_mf.txt", GMT_MF)
np.savetxt("GMT/gmt_wf.txt", GMT_WF)
np.savetxt("GMT/gmt_sf.txt", GMT_SF)





fig, (ax0, ax1) = plt.subplots(nrows=1, ncols=2, figsize=(16,6))

ax0.plot(np.arange(0, 1000, 1), x_SF[:1000], color="red", label="Business as usual")
ax0.plot(np.arange(0, 1000, 1), x_MF[:1000], color="blue", label="Moderately incentivized clean energy")
ax0.plot(np.arange(0, 1000, 1), x_WF[:1000], color="green", label="Highly incentivized clean energy")
ax0.set_xlabel("Time [yrs]")
ax0.set_ylabel("CO2 [ppm]")
ax0.set_ylim([300, 900])
ax0.legend(loc="upper right", fontsize=15)


ax1.plot(np.arange(0, 1000, 1), GMT_SF[:1000], color="red", label="Business as usual")
ax1.plot(np.arange(0, 1000, 1), GMT_MF[:1000], color="blue", label="Moderately incentivized clean energy")
ax1.plot(np.arange(0, 1000, 1), GMT_WF[:1000], color="green", label="Highly incentivized clean energy")
ax1.set_xlabel("Time [yrs]")
ax1.set_ylabel("$\Delta$ GMT [°C]")
ax1.set_ylim([0.0, 3.0])
ax1.legend(loc="upper right", fontsize=15)

fig.tight_layout()
fig.savefig("co2_gmt.png")
fig.savefig("co2_gmt.pdf")
# fig.show()
fig.clf()
plt.close()


print("Finish!")