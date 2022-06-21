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





#CO2 values
x_MF = np.arange(280, 1000, 1)

#Conversion to GMT
GMT_MF = []
for i in range(0, len(x_MF)):
    gmt_val_MF = co2_gmt(x_MF[i])
    GMT_MF.append(gmt_val_MF)
print(GMT_MF[0])
GMT_MF = np.array(GMT_MF)/2.0226354355276586 #the 2.0226354355276586-factor is the rescaling factor to set 400ppm to 1.0°C above pre-industrial

for i in range(0, len(x_MF)):
    print(x_MF[i], GMT_MF[i])


fig, (ax0, ax1, ax2) = plt.subplots(nrows=1, ncols=3, figsize=(24,6))

ax0.plot(np.arange(0, len(x_MF), 1), x_MF, color="red")
ax0.set_xlabel("Time [yrs]")
ax0.set_ylabel("CO2 [ppm]")
#ax0.set_ylim([300, 900])

ax1.plot(np.arange(0, len(x_MF), 1), GMT_MF, color="red")
ax1.set_xlabel("Time [yrs]")
ax1.set_ylabel("$\Delta$ GMT [°C]")
#ax1.set_ylim([0.0, 3.0])

ax2.plot(x_MF, GMT_MF, color="blue")
ax2.set_xlabel("CO [ppm]")
ax2.set_ylabel("$\Delta$ GMT [°C]")

fig.tight_layout()
fig.savefig("co2_gmt_test.png")
fig.savefig("co2_gmt_test.pdf")
# fig.show()
fig.clf()
plt.close()


print("Finish!")