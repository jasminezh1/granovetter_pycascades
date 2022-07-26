# Add modules directory to path
from audioop import avg
import os
import sys

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
duration = 6000 
# previously 50000


#Names to create the respective directories
namefile = "no"
if switch_unc_barrett == False and switch_unc == True:
    long_save_name = "results_tempRateChanges"
elif switch_unc_barrett == True and switch_unc == True:
    long_save_name = "results_barrett_3"
elif switch_unc_barrett == False and switch_unc == False:
    long_save_name = "results_no_uncertainty_3"
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
a = 0.14
ets = 0.1
c = 0.6


################################# MAIN LOOP #################################
#for kk in plus_minus_links:
finalTemp = []
finalTipped = []
finalActive = []
finalTimes = []

kk = [-1,0,-1]
#EnvToSocCoupling = np.linspace(0.0, 1, 21, endpoint=True)
tempRate = np.linspace(0.005, 0.05, 21, endpoint=True)

for tr in tempRate:
    allTipped = False

    #print("Wais to Thc:{}".format(kk[0]))
    #print("Amaz to Nino:{}".format(kk[1]))
    #print("Thc to Amaz:{}".format(kk[2]))
    try:
        os.stat("{}".format(long_save_name))
    except:
        os.makedirs("{}".format(long_save_name))

    try:
        os.stat("{}/{}_feedbacks".format(long_save_name, namefile))
    except:
        os.mkdir("{}/{}_feedbacks".format(long_save_name, namefile))

    try:
        os.stat("{}/{}_feedbacks/{:.2f}_{:.2f}".format(long_save_name, namefile, a, c))
    except:
        os.mkdir("{}/{}_feedbacks/{:.2f}_{:.2f}".format(long_save_name, namefile, a, c))

    try:
        os.stat("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))
    except:
        os.mkdir("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))

    try:
        os.stat("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))
    except:
        os.mkdir("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))

    #save starting conditions
    np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/empirical_values_care{}_{}.txt".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2], a, c), sys_var, delimiter=" ", fmt="%s")


    
    times = []

    #for strength in coupling_strength:
    strength = 0.2
    #print("Coupling strength: {:.2f}".format(strength))

    print("tr: ", tr)

        
    # if os.path.isfile("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_{:.3f}.txt".format(long_save_name, 
    # namefile, ets, c, kk[0], kk[1], kk[2], ets, c, strength, ets)) == True:
    #     print("File already computed")
    #     continue

    tippedElements = 0
    granTipped = False
    output = []
    gmt = [0]
    # roots = np.empty((duration,3,))
    # roots[:] = np.nan
    first = True
    firstRoot = []
    numTipped = []
    closeGis, closeThc, closeWais, closeAmaz = -1.0, -1.0, -1.0, -1.0

    for t in range(2, int(duration)+2):
        

        ######THE SOCIAL MODEL
        #model = gwmodel(threshold_frac,avg_degree,a+tippedElements*0.05,c-tippedElements*0.05)

        # take the negative ones (not yet tipped) between -1 and 0
        # only care about specific value for negatives (cuz that's how close)
        # the smaller the number the better
        changeGis = np.min(closeGis, 0)
        changeThc = np.min(closeThc, 0)
        changeWais = np.min(closeWais, 0)
        changeAmaz = np.min(closeAmaz, 0)
        totalChange = changeGis + changeThc + changeWais + changeAmaz + 4  # total change is positive, 0 to 4
        totalChange = totalChange / 4 * c * ets

        model = gwmodel(threshold_frac,avg_degree, a+totalChange, c-totalChange)
        ###END SOCIAL MODEL

        ######THE NATURAL MODEL
        effective_GMT = gmt[-1]

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
        closeGis = ev.get_timeseries()[1][-1, 0]
        closeThc = ev.get_timeseries()[1][-1, 1], 
        closeWais = ev.get_timeseries()[1][-1, 2]
        closeAmaz = ev.get_timeseries()[1][-1, 3]

        if((not allTipped) and tippedElements==4):
            finalTimes.append(t-1)
            allTipped = True

        if t==2:
            root_x, root_y = model.guess(a)
        else:
            root_x, root_y = model.guess(firstRoot[-1])

        activeShare = float(root_x[0])
        firstRoot.append(activeShare)
        numTipped.append(tippedElements / 4)

        lastTemp = float(gmt[-1])
        # stop the temperature from increasing forever
        if(activeShare > (a+c-0.05)):
            activeShare = 1             # in the future can make >1
            #newTemp = float(lastTemp - 0.0005)
        #else:
        newTemp = float(lastTemp + (1 - activeShare) * tr)        # vary this 0.01 from 0.005 to 0.05
        #print(newTemp," ",type(newTemp))
        gmt.append(newTemp)

        if(len(root_x)==1):
            granTipped = True
            #  unclear but sometimes list sometimes scalar
            if(isinstance(root_x, list)):
                root_x = root_x[0]
            #roots[t-2][0] = float(root_x)       # used to be root_x[0]
            if(first):
                times.append(t-2)
            first = False
        #else:
            #roots[t-2] = root_x

        #saving structure
        output.append([t,
                        ev.get_timeseries()[1][-1, 0],
                        ev.get_timeseries()[1][-1, 1],
                        ev.get_timeseries()[1][-1, 2],
                        ev.get_timeseries()[1][-1, 3],
                        net.get_number_tipped(ev.get_timeseries()[1][-1]),
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

        saveGMT = data[-1]
        if(isinstance(saveGMT[-2],np.ndarray)):
            saveGMT[-2] = saveGMT[-2][0]

        np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.3f}.txt".format(long_save_name, namefile, a, c, 
            kk[0], kk[1], kk[2], a, c, strength, ets, tr), saveGMT)
        time = data.T[0]
        state_gis = data.T[1]
        state_thc = data.T[2]
        state_wais = data.T[3]
        state_amaz = data.T[4]
        root_series = data.T[-2] #roots ????
        gmt_series = data.T[-1] #gmt series

        np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.3f}_gmt.txt".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength, ets, tr), gmt_series)
        
        # np.savetxt("{}/{}_feedbacks/{}_{}/network_{}_{}_{}/feedbacks_care{}_{}_{:.2f}_{:.2f}_ALLROOTS.txt".format(long_save_name, namefile, a, c,
        #         kk[0], kk[1], kk[2], a, c, strength, ets), roots)

        #plotting structure
        fig = plt.figure()
        plt.grid(True)
        plt.title("Temp Rate: {:.3f} A: {:.2f} C: {:.2f} ETS: {:.2f}\n  Wais to Thc:{}  Amaz to Nino:{} Thc to Amaz:{}".format(
            tr, a, c, ets, kk[0], kk[1], kk[2]))
        plt.plot(time, state_gis, label="GIS", color='c')
        plt.plot(time, state_thc, label="THC", color='b')
        plt.plot(time, state_wais, label="WAIS", color='k')
        plt.plot(time, state_amaz, label="AMAZ", color='g')
        plt.xlabel("Time [yr]")
        plt.ylabel("system feature f [a.u.]")
        plt.legend(loc='best')  # , ncol=5)
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.3f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength, ets, tr))
        #plt.show()
        plt.clf()
        plt.close()


        fig = plt.figure()
        plt.grid(True)
        plt.title("Temp Rate: {:.3f} A: {:.2f} C: {:.2f} ETS: {:.2f}\n  Wais to Thc:{}  Amaz to Nino:{} Thc to Amaz:{}".format(
            tr, a, c, ets, kk[0], kk[1], kk[2]))
        plt.plot(time, gmt_series, label="GMT", color='r')
        plt.xlabel("Time [yr]")
        plt.ylabel("$\Delta$ GMT [°C]")
        plt.legend(loc='best')  # , ncol=5)
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/gmtseries_care{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.3f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength, ets, tr))
        #plt.show()
        plt.clf()
        plt.close()

        
        # **********************************
        fig, ax1 = plt.subplots()
        plt.grid(True)
        plt.title("Temp Rate {:.3f} A: {:.2f} C: {:.2f} ETS: {:.2f}\n  Wa to Thc:{}  Am to Ni:{} Thc to Am:{}".format(
            tr, a, c, ets, kk[0], kk[1], kk[2]))
        ax1.set_xlabel("Time [yr]")
        ax1.set_ylabel("Percentage")
        plot1 = ax1.plot(time, firstRoot, label = "active people", color='r')
        plot2 = ax1.plot(time, numTipped, label = "tipped elements", color='g')

        ax2 = ax1.twinx()
        ax2.set_ylabel("GMT", color = 'b')
        plot3 = ax2.plot(time, gmt_series, label="GMT", color='b')
        #plt.legend(loc='best')
        lns = plot1 + plot2 + plot3
        labels = [l.get_label() for l in lns]
        plt.legend(lns, labels, loc=0)
        ax1.set_ylim(top = 1.25)
        ax2.set_ylim(top = 5)
        fig.tight_layout()
        plt.savefig("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/root_care{:.2f}_{:.2f}_{:.2f}_{:.2f}_{:.3f}.pdf".format(long_save_name, namefile, a, c,
            kk[0], kk[1], kk[2], a, c, strength, ets, tr))
        #plt.show()
        plt.clf()
        plt.close()

        
        
    #print("times: ", times)



    finalTemp.append(gmt_series[-1])
    finalTipped.append(numTipped[-1])
    finalActive.append(firstRoot[-1])

#it is necessary to limit the amount of saved files
#--> compose one pdf file for each network setting and remove the other time-files
current_dir = os.getcwd()
os.chdir("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/".format(long_save_name, namefile, a, c, kk[0], kk[1], kk[2]))
pdfs = np.array(np.sort(glob.glob("feedbacks_care{:.2f}_{:.2f}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm feedbacks_care{:.2f}_{:.2f}_*.pdf".format(a, c))
    merger.write("feedbacks_complete_care{:.2f}_{:.2f}.pdf".format(a, c))
    print("Complete PDFs merged - Part 1")

pdfs = np.array(np.sort(glob.glob("gmtseries_care{:.2f}_{:.2f}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm gmtseries_care{:.2f}_{:.2f}_*.pdf".format(a,c))
    merger.write("gmtseries_complete_care{:.2f}_{:.2f}.pdf".format(a, c))
    print("Complete PDFs merged - Part 2")

pdfs = np.array(np.sort(glob.glob("root_care{:.2f}_{:.2f}_*.pdf".format(a, c)), axis=0))
if len(pdfs) != 0.:
    merger = PdfFileMerger()
    for pdf in pdfs:
        merger.append(pdf)
    os.system("rm root_care{:.2f}_{:.2f}_*.pdf".format(a, c))
    merger.write("root_complete_care{:.2f}_{:.2f}.pdf".format(a, c))
    print("Complete PDFs merged - Part 3")

os.chdir(current_dir)


np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_finalTemps.txt".format(long_save_name, namefile, a, c,
    kk[0], kk[1], kk[2], a, c, strength), finalTemp)  
np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_finalTipped.txt".format(long_save_name, namefile, a, c,
    kk[0], kk[1], kk[2], a, c, strength), finalTipped)  
np.savetxt("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/feedbacks_care{:.2f}_{:.2f}_{:.2f}_finalActive.txt".format(long_save_name, namefile, a, c,
    kk[0], kk[1], kk[2], a, c, strength), finalActive) 

fig, ax1 = plt.subplots()
plt.grid(True)
plt.title("A: {:.2f} C: {:.2f} ETS: {:.2f}\n  Wa to Thc:{}  Am to Ni:{} Thc to Am:{}".format(
    a, c, ets, kk[0], kk[1], kk[2]))
ax1.set_xlabel("Rate of Temperature Change from Active People")
ax1.set_ylabel("Time All Elements Tipped")
plot1 = ax1.plot(tempRate, finalTimes, label = "time", color='r')
#plot2 = ax1.plot(tempRate, finalTipped, label = "tipped elements", color='g')

ax2 = ax1.twinx()
ax2.set_ylabel("GMT", color = 'b')
plot3 = ax2.plot(tempRate, finalTemp, label="GMT", color='b')
#plt.legend(loc='best')
lns = plot1 + plot3
labels = [l.get_label() for l in lns]
plt.legend(lns, labels, loc=0)
#ax1.set_ylim(top = 1.25)
#ax2.set_ylim(top = 5)
fig.tight_layout()
plt.savefig("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/ALL2_{:.2f}_{:.2f}_{:.2f}.pdf".format(long_save_name, namefile, a, c,
    kk[0], kk[1], kk[2], a, c, strength))
#plt.show()
plt.clf()
plt.close()   

fig, ax1 = plt.subplots()
plt.grid(True)
plt.title("A: {:.2f} C: {:.2f} ETS: {:.2f}\n  Wa to Thc:{}  Am to Ni:{} Thc to Am:{}".format(
    a, c, ets, kk[0], kk[1], kk[2]))
ax1.set_xlabel("Rate of Temperature Change from Active People")
ax1.set_ylabel("Percentage")
plot1 = ax1.plot(tempRate, finalActive, label = "active people", color='r')
plot2 = ax1.plot(tempRate, finalTipped, label = "tipped elements", color='g')

ax2 = ax1.twinx()
ax2.set_ylabel("GMT", color = 'b')
plot3 = ax2.plot(tempRate, finalTemp, label="GMT", color='b')
#plt.legend(loc='best')
lns = plot1 + plot2 + plot3
labels = [l.get_label() for l in lns]
plt.legend(lns, labels, loc=0)
#ax1.set_ylim(top = 1.25)
#ax2.set_ylim(top = 5)
fig.tight_layout()
plt.savefig("{}/{}_feedbacks/{:.2f}_{:.2f}/network_{}_{}_{}/ALL_{:.2f}_{:.2f}_{:.2f}.pdf".format(long_save_name, namefile, a, c,
    kk[0], kk[1], kk[2], a, c, strength))
#plt.show()
plt.clf()
plt.close()   

print("Finish")
#end = time.time()
#print("Time elapsed until Finish: {}s".format(end - start))
