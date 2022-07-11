# Add modules directory to path
from audioop import avg
import os
import sys
import re

sys.path.append('')

# global imports
import numpy as np
import matplotlib
#matplotlib.use("Agg")
matplotlib.use( 'tkagg' )
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.)
import itertools
import time
import glob
from PyPDF2 import PdfFileMerger
from netCDF4 import Dataset
import tkinter

# private imports from sys.path
from core.evolve import evolve

#private imports for earth system
from earth_sys.timing_no_enso import timing
from earth_sys.functions_earth_system_no_enso import global_functions
from earth_sys.earth_no_enso import earth_system

#private imports from Saverio's model
from earth_sys.nzmodel import nzmodel
from gran_wied.gwmodel import gwmodel


#for cluster computations
os.chdir("/Users/jasminezhang/Documents/PIK22")

#measure time
#start = time.time()
#############################GLOBAL SWITCHES#########################################
time_scale = True            # time scale of tipping is incorporated
plus_minus_include = True    # from Kriegler, 2009: Unclear links; if False all unclear links are set to off state and only network "0-0" is computed
switch_unc = True            # if uncertainty is taken into account at all
switch_unc_barrett = False   # is there an uncertainty in the thresholds, which lowers the ability to cooperate (following Barrett et al., Nature Climate Change, 2014)
######################################################################


# switch_unc_barrett cannot be True in case switch_unc is False
if switch_unc == False:
    if switch_unc_barrett == True:
        print("Switch_unc_barrett cannot be True in case switch_unc is False - Please check !!!")
        print("The percevied uncertainty (Barrett) does not make sense if overall uncertainty about threshold is zero")
        quit()

## temperature uncertainty that has to be taken into account
unc_temp = np.array([[0.8, 3.0], [1.4, 8.0], [1.0, 3.0], [2.0, 6.0]]) #limit of critical temperature thresholds of GIS, THC, WAIS, AMAZ


#actual real simulation years
duration = 1000 
# previously 50000


#Names to create the respective directories
namefile = "no"
if switch_unc_barrett == False and switch_unc == True:
    long_save_name = "results_2"
elif switch_unc_barrett == True and switch_unc == True:
    long_save_name = "results_barrett_2"
elif switch_unc_barrett == False and switch_unc == False:
    long_save_name = "results_no_uncertainty_2"
    print("We don't want to compute this at the moment!")
    quit()
else:
    print("Check switches of uncertainty")
    quit()

#######################GLOBAL VARIABLES##############################
#drive coupling strength
coupling_strength = np.linspace(0.0, 1.0, 6, endpoint=True)
# previusly 11

########################Declaration of variables from passed values#######################
#Must sort out first value since this is the actual file 
sys_var = np.array(sys.argv[1:], dtype=str) #low sample -3, intermediate sample: -2, high sample: -1
#####################################################################


#Tipping ranges from distribution
limits_gis, limits_thc, limits_wais, limits_amaz, limits_nino = float(sys_var[0]), float(sys_var[1]), float(sys_var[2]), float(sys_var[3]), float(sys_var[4])

#Probability fractions
# TO GIS
pf_wais_to_gis, pf_thc_to_gis = float(sys_var[5]), float(sys_var[6])
# TO THC
pf_gis_to_thc, pf_nino_to_thc, pf_wais_to_thc = float(sys_var[7]), float(sys_var[8]), float(sys_var[9])
# TO WAIS
pf_nino_to_wais, pf_thc_to_wais, pf_gis_to_wais = float(sys_var[10]), float(sys_var[11]), float(sys_var[12])
# TO NINO
pf_thc_to_nino, pf_amaz_to_nino = float(sys_var[13]), float(sys_var[14])
# TO AMAZ
pf_nino_to_amaz, pf_thc_to_amaz = float(sys_var[15]), float(sys_var[16])

#tipping time scales
tau_gis, tau_thc, tau_wais, tau_nino, tau_amaz = float(sys_var[17]), float(sys_var[18]), float(sys_var[19]), float(sys_var[20]), float(sys_var[21])



#Time scale
"""
All tipping times are computed ion comparison to the Amazon rainforest tipping time. As this is variable now, this affects the results to a (very) level
"""
if time_scale == True:
    print("compute calibration timescale")
    #function call for absolute timing and time conversion
    time_props = timing(tau_gis, tau_thc, tau_wais, tau_amaz, tau_nino)
    gis_time, thc_time, wais_time, nino_time, amaz_time = time_props.timescales()
    conv_fac_gis = time_props.conversion()
else:
    #no time scales included
    gis_time = thc_time = wais_time = nino_time = amaz_time = 1.0
    conv_fac_gis = 1.0

#include uncertain "+-" links:
if plus_minus_include == True:
    plus_minus_links = np.array(list(itertools.product([-1.0, 0.0, 1.0], repeat=3)))

    #in the NO_ENSO case (i.e., the second link must be 0.0)
    plus_minus_data = []
    for pm in plus_minus_links:
        if pm[1] == 0.0:
            plus_minus_data.append(pm)
    plus_minus_links = np.array(plus_minus_data)

else:
    plus_minus_links = [np.array([1., 1., 1.])]


#directories for the Monte Carlo simulation
#mc_dir = int(sys_var[-1])
mc_dir = 0


################################# MAIN #################################
#Create Earth System
earth_system = earth_system(gis_time, thc_time, wais_time, nino_time, amaz_time,
                            limits_gis, limits_thc, limits_wais, limits_nino, limits_amaz,
                            pf_wais_to_gis, pf_thc_to_gis, pf_gis_to_thc, pf_nino_to_thc,
                            pf_wais_to_thc, pf_gis_to_wais, pf_thc_to_wais, pf_nino_to_wais,
                            pf_thc_to_nino, pf_amaz_to_nino, pf_nino_to_amaz, pf_thc_to_amaz)

#####Variables to be set
#Create Social Model

threshold_frac = 0.5
avg_degree = 10
a = 0.23
c = 0.5


################################# MAIN LOOP #################################
#for kk in plus_minus_links:
kk = [-1,0,-1]

print("Wais to Thc:{}".format(kk[0]))
print("Amaz to Nino:{}".format(kk[1]))
print("Thc to Amaz:{}".format(kk[2]))
try:
    os.stat("{}".format(long_save_name))
except:
    os.makedirs("{}".format(long_save_name))

try:
    os.stat("{}/{}_feedbacks".format(long_save_name, namefile))
except:
    os.mkdir("{}/{}_feedbacks".format(long_save_name, namefile))

try:
    os.stat("{}/{}_feedbacks/{}_{}".format(long_save_name, namefile, a, c))
except:
    os.mkdir("{}/{}_feedbacks/{}_{}".format(long_save_name, namefile, a, c))

try:
    os.stat("{}/{}_feedbacks/{}_{}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))
except:
    os.mkdir("{}/{}_feedbacks/{}_{}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))


#save starting conditions
np.savetxt("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/empirical_values_care{}_{}.txt".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2], a, c), sys_var, delimiter=" ", fmt="%s")

times = []

for strength in coupling_strength:
    print("Coupling strength: {:.2f}".format(strength))
    tippedElements = 0
    granTipped = False
    output = []
    gmt = [0]
    roots = np.empty((duration,3,))
    roots[:] = np.nan
    first = True
    firstRoot = []
    numTipped = []

    for t in range(2, int(duration)+2):
        
        if os.path.isfile("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}.txt".format(long_save_name, 
            namefile, a, c, kk[0], kk[1], kk[2], a, c, strength)) == True:
            print("File already computed")
            break
        
        ######THE SOCIAL MODEL
        model = gwmodel(threshold_frac,avg_degree,a+tippedElements*0.05,c-tippedElements*0.05)
        ###END SOCIAL MODEL

        ######THE NATURAL MODEL
        effective_GMT = gmt[-1]
        #print(effective_GMT)

        #get back the network of the Earth system
        net = earth_system.earth_network(effective_GMT, strength, kk[0], kk[1], kk[2])

        # initialize state
        if t == 2:
            initial_state = [-1, -1, -1, -1] #initial state
        else:
            initial_state = [ev.get_timeseries()[1][-1, 0], ev.get_timeseries()[1][-1, 1], ev.get_timeseries()[1][-1, 2], ev.get_timeseries()[1][-1, 3]]
        ev = evolve(net, initial_state)
        # plotter.network(net)

        # Timestep to integration; it is also possible to run integration until equilibrium
        timestep = 0.1

        #t_end given in years; also possible to use equilibrate method
        t_end = 1.0/conv_fac_gis #simulation length in "real" years, in this case 1 year
        ev.integrate(timestep, t_end)
        #######END: THE NATURAL MODEL

        tippedElements = net.get_number_tipped(ev.get_timeseries()[1][-1])
        root_x, root_y = model.guess()

        activeShare = float(root_x[0])
        firstRoot.append(activeShare)
        numTipped.append(tippedElements / 4)

        # stop the temperature from increasing forever
        if(activeShare > (a+c-0.05)):                   # ********** prematurely stops temp growth
            activeShare = 1

        lastTemp = float(gmt[-1])
        newTemp = float(lastTemp + (1 - activeShare) * 0.05)
        gmt.append(newTemp)

        if(len(root_x)==1):
            granTipped = True
            #  unclear but sometimes list sometimes scalar
            if(isinstance(root_x, list)):
                root_x = root_x[0]
            roots[t-2][0] = float(root_x)       # used to be root_x[0]
            if(first):
                times.append(t-2)
                #print(root_x)
            first = False
        else:
            roots[t-2] = root_x

        #saving structure
        output.append([t,
                        ev.get_timeseries()[1][-1, 0],
                        ev.get_timeseries()[1][-1, 1],
                        ev.get_timeseries()[1][-1, 2],
                        ev.get_timeseries()[1][-1, 3],
                        net.get_number_tipped(ev.get_timeseries()[1][-1]),
                        # last state [0 - 3]. states are an array? maybe state of each element?
                        [net.get_tip_states(ev.get_timeseries()[1][-1])[0]].count(True),   
                        [net.get_tip_states(ev.get_timeseries()[1][-1])[1]].count(True),
                        [net.get_tip_states(ev.get_timeseries()[1][-1])[2]].count(True),
                        [net.get_tip_states(ev.get_timeseries()[1][-1])[3]].count(True),
                        root_x, effective_GMT
                        ])

    
    #necessary for break condition
    if len(output) != 0:
        #saving structure
        data = np.array(output, dtype=object)

        np.savetxt("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}.txt".format(long_save_name, namefile, a, c, 
            kk[0], kk[1], kk[2], a, c, strength), data[-1])
        time = data.T[0]
        state_gis = data.T[1]
        state_thc = data.T[2]
        state_wais = data.T[3]
        state_amaz = data.T[4]
        root_series = data.T[-2] #roots ????
        gmt_series = data.T[-1] #gmt series

        np.savetxt("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}_gmt.txt".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength), gmt_series)
        
        np.savetxt("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}_ALLROOTS.txt".format(long_save_name, namefile, a, c,
                kk[0], kk[1], kk[2], a, c, strength), roots)

        #plotting structure
        fig = plt.figure()
        plt.grid(True)
        plt.title("Coupling strength: {}\n  Wais to Thc:{}  Amaz to Nino:{} Thc to Amaz:{}".format(
            np.round(strength, 2), kk[0], kk[1], kk[2]))
        plt.plot(time, state_gis, label="GIS", color='c')
        plt.plot(time, state_thc, label="THC", color='b')
        plt.plot(time, state_wais, label="WAIS", color='k')
        plt.plot(time, state_amaz, label="AMAZ", color='g')
        plt.xlabel("Time [yr]")
        plt.ylabel("system feature f [a.u.]")
        plt.legend(loc='best')  # , ncol=5)
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength))
        #plt.show()
        plt.clf()
        plt.close()


        fig = plt.figure()
        plt.grid(True)
        plt.title("Coupling strength: {}\n  Wais to Thc:{}  Amaz to Nino:{} Thc to Amaz:{}".format(
            np.round(strength, 2), kk[0], kk[1], kk[2]))
        plt.plot(time[0:2000], gmt_series[0:2000], label="GMT", color='r')
        plt.xlabel("Time [yr]")
        plt.ylabel("$\Delta$ GMT [°C]")
        plt.legend(loc='best')  # , ncol=5)
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/gmtseries_care{}_{}_{:.2f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength))
        #plt.show()
        plt.clf()
        plt.close()

        
        # **********************************
        fig = plt.figure()
        plt.grid(True)
        plt.title("Coupling strength: {}  A: {}  C: {}\n  Wa to Thc:{}  Am to Ni:{} Thc to Am:{}".format(
            np.round(strength, 2), a, c, kk[0], kk[1], kk[2]))
        #plt.title("Active People and Tipped Elements \nfor initial a = {} c = {}".format(a,c))
        plt.plot(time[0:duration], firstRoot[0:duration], label = "active people", color='r')
        plt.plot(time[0:duration], numTipped[0:duration], label = "tipped elements", color='g')
        plt.xlabel("Time [yr]")
        plt.ylabel("Percentage")
        plt.legend(loc='best')
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/root_care{}_{}_{:.2f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength))
        #plt.show()
        plt.clf()
        plt.close()

    
    
#print("times: ", times)

# it is necessary to limit the amount of saved files
# --> compose one pdf file for each network setting and remove the other time-files
current_dir = os.getcwd()
os.chdir("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))
pdfs = np.array(np.sort(glob.glob("feedbacks_care{}_{}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm feedbacks_care{}_{}_*.pdf".format(a, c))
    merger.write("feedbacks_complete_care{}_{}.pdf".format(a, c))
    print("Complete PDFs merged - Part 1")

pdfs = np.array(np.sort(glob.glob("gmtseries_care{}_{}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm gmtseries_care{}_{}_*.pdf".format(a, c))
    merger.write("gmtseries_complete_care{}_{}.pdf".format(a, c))
    print("Complete PDFs merged - Part 2")

pdfs = np.array(np.sort(glob.glob("root_care{}_{}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm root_care{}_{}_*.pdf".format(a, c))
    merger.write("root_complete_care{}_{}.pdf".format(a, c))
    print("Complete PDFs merged - Part 3")

os.chdir(current_dir)

    

print("Finish")
#end = time.time()
#print("Time elapsed until Finish: {}s".format(end - start))
