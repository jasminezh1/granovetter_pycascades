# Add modules directory to path
import os
import sys
import re

sys.path.append('')

# global imports
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
sns.set(font_scale=1.)
import itertools
import time
import glob
from PyPDF2 import PdfFileMerger
from netCDF4 import Dataset


# private imports from sys.path
from core.evolve import evolve

#private imports for earth system
from earth_sys.timing_no_enso import timing
from earth_sys.functions_earth_system_no_enso import global_functions
from earth_sys.earth_no_enso import earth_system

#private imports from Saverio's model
from earth_sys.nzmodel import nzmodel



#for cluster computations
#os.chdir("/p/projects/dominoes/nicowun/conceptual_tipping/uniform_distribution/socio_climate")
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
duration = 50000 



#Names to create the respective directories
namefile = "no"
if switch_unc_barrett == False and switch_unc == True:
    long_save_name = "results"
elif switch_unc_barrett == True and switch_unc == True:
    long_save_name = "results_barrett"
elif switch_unc_barrett == False and switch_unc == False:
    long_save_name = "results_no_uncertainty"
    print("We don't want to compute this at the moment!")
    quit()
else:
    print("Check switches of uncertainty")
    quit()

#######################GLOBAL VARIABLES##############################
#drive coupling strength
coupling_strength = np.linspace(0.0, 1.0, 11, endpoint=True)


"""
#make ensemble a bit smaller for the moment
GMT_files_1 = np.sort(glob.glob("temp_input/timeseries_final/*100.txt"))
GMT_files_2 = np.sort(glob.glob("temp_input/timeseries_final/*500.txt"))
GMT_files_3 = np.sort(glob.glob("temp_input/timeseries_final/*1000.txt"))
GMT_files = np.concatenate((np.concatenate((GMT_files_1, GMT_files_2)), GMT_files_3))
"""


########################Declaration of variables from passed values#######################
#Must sort out first and second value since this is the actual file and the number of nodes used
sys_var = np.array(sys.argv[1:], dtype=str) #low sample -3, intermediate sample: -2, high sample: -1
#print("start ", sys.argv[1])

#TEST SYTEM
#print("USING TEST SYSTEM")
#sys_var = np.array([1.6, 4.75, 3.25, 4.0, 4.0, 0.2, 1.0, 1.0, 0.2, 0.3, 0.5, 0.15, 1.0, 0.2, 0.15, 1.0, 0.4, 4000, 150, 4000, 50, 100, 2400])
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

#care (in fraction of mu*E_CL)
care = float(sys_var[22])
#tau
tau = float(sys_var[23]) #tau_low=tau/2.2; tau_high=tau*2.2; get low, indermediate, high scenario
if tau > 15 and tau < 17:
    scenario = "sf"
elif tau > 35 and tau < 39:
    scenario = "mf"
elif tau > 80 and tau < 82:
    scenario = "wf"
else:
    print("Outside boundaries - wrong tau value - please CHECK!!!")
    quit()



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
mc_dir = int(sys_var[-1])

#plus_minus_links = [plus_minus_links[0]]
#plus_minus_links = [plus_minus_links[1]]
#plus_minus_links = [plus_minus_links[2]]
#plus_minus_links = [plus_minus_links[3]]
#plus_minus_links = [plus_minus_links[4]]
#plus_minus_links = [plus_minus_links[5]]
#plus_minus_links = [plus_minus_links[6]]
#plus_minus_links = [plus_minus_links[7]]
#plus_minus_links = [plus_minus_links[8]]



#plus_minus_links = plus_minus_links[6:8]
#print(plus_minus_links)


################################# MAIN #################################
#Create Earth System
earth_system = earth_system(gis_time, thc_time, wais_time, nino_time, amaz_time,
                            limits_gis, limits_thc, limits_wais, limits_nino, limits_amaz,
                            pf_wais_to_gis, pf_thc_to_gis, pf_gis_to_thc, pf_nino_to_thc,
                            pf_wais_to_thc, pf_gis_to_wais, pf_thc_to_wais, pf_nino_to_wais,
                            pf_thc_to_nino, pf_amaz_to_nino, pf_nino_to_amaz, pf_thc_to_amaz)

#####Variables to be set
#Create Social Model
NS=12.14 #(* ppm reduction in CO2 conc. due to natural sinks for the year 2018*)
a=1
OF=0
E0=171.24 #(*10^6 GW Total energy demand for the year 2018*)
D = 0.0
# Eta=1.23688155922039 #(* ppm/GW increase in CO2 conc. per GW of energy produced from fossil fuels*)*)
Kappa1=0.08
Kappa2=0.2 #(*2.25*)
K2=0.2
# (*K3=0.05;*)
K3=0.2
n=4 #number of blocks
yin=17.77 # (*10^6 GW Clean energy production for the year 2018*)
zin=330 #(*$billion Clean energy Investment for the year 2018*)
K1=(0.1759+K2*yin)/zin
xin=280.0    #408.52
Eta=(2.3+NS+OF)/(E0-yin) # (* ppm/GW increase in CO2 conc. per GW of energy produced from fossil fuels*)

a1=80
a2=0.05
a3=18.2
a4=0.0027
xref=280 #ppm pre-industrial CO2 concentration
ScF=1
phi=1.4231*ScF
phi1=phi
phi2=phi
phi3=phi
phi4=phi
epsilon=np.multiply(0.236*1.4231, 0.236*1.4231)*ScF
epsilon1=epsilon
epsilon2=epsilon
epsilon3=epsilon
epsilon4=epsilon
# tau1=30
tau1=tau #37
tau11=tau1
tau12=tau1
tau13=tau1
tau14=tau1
# mu=.005
mu=.0035
mu1=mu
mu2=mu
mu3=mu
mu4=mu





nzmodel = nzmodel(NS, a, OF, Kappa1, Kappa2, K2, K3, n, K1, 
        Eta, a1, a2, a3, a4, xref, ScF, phi, phi1, phi2, phi3, phi4, 
        epsilon, epsilon1, epsilon2, epsilon3, epsilon4, 
        tau1, tau11, tau12, tau13, tau14, mu, mu1, mu2, mu3, mu4)



################################# MAIN LOOP #################################
for kk in plus_minus_links:
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
        os.stat("{}/{}_feedbacks/network_{}_{}_{}".format(long_save_name, namefile, kk[0], kk[1], kk[2]))
    except:
        os.mkdir("{}/{}_feedbacks/network_{}_{}_{}".format(long_save_name, namefile, kk[0], kk[1], kk[2]))

    try:
        os.stat("{}/{}_feedbacks/network_{}_{}_{}/{}".format(long_save_name, namefile, kk[0], kk[1], kk[2], str(mc_dir).zfill(4) ))
    except:
        os.mkdir("{}/{}_feedbacks/network_{}_{}_{}/{}".format(long_save_name, namefile, kk[0], kk[1], kk[2], str(mc_dir).zfill(4) ))

    #save starting conditions
    np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/empirical_values_care{:.2f}_{}.txt".format(long_save_name, namefile, kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario), sys_var, delimiter=" ", fmt="%s")

    for strength in coupling_strength:
        print("Coupling strength: {}".format(strength))

        ###DEFINE START VALUES FOR SOCIAL MODEL CORRECTLY
        #reference values to be used in the following FOR-LOOP
        E0_ref = 171.24
        E0 = 171.24         #E: Energy demand at t=0
        xin=280.0           #dCO2/dt: CO2 concentration at t=0 [Note that we start at preindustrial levels of CO2 and not at 408 ppm as in 2018]
        yin=17.77           #E_CL: clean energy production at t=0
        zin=330             #SA: socio-political acceptability
        D = 0.0             #PED: perceived economic damage PED
        effective_GMT = 0   #GMT: Temperature at t=0
        ####


        output = []
        for t in range(2, int(duration)+2):
            #print(t)
            if os.path.isfile("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}.txt".format(long_save_name, 
                namefile, kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength)) == True:
                print("File already computed")
                break
            


            ######THE SOCIAL MODEL
            limit_closest = np.amin([limits_gis, limits_thc, limits_wais, limits_amaz]) #closest critical threshold
            ind_closest = np.argmin([limits_gis, limits_thc, limits_wais, limits_amaz]) #index of closest critical threshold to later identify the tipping element
            gmt, E_arr, x_arr, y_arr, z_arr, D_arr = nzmodel.NetZeroFunDamage(t, E0, E0_ref, xin, yin, zin, D, effective_GMT, limit_closest, ind_closest, switch_unc, switch_unc_barrett, care, unc_temp)
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


            E0 = E_arr[-1]
            xin = x_arr[-1]
            yin = y_arr[-1]
            zin = z_arr[-1]
            D = D_arr[-1]

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
                           E0, xin, yin, zin, D, effective_GMT
                           ])



        """
        output_test = np.array(output)
        gmt = output_test.T[-1]
        co2 = output_test.T[-5]
        print(co2)
        print(gmt)
        
        #Plotting
        fig, ((ax0)) = plt.subplots(nrows=1, ncols=1, figsize=(10, 8))

        T_final = 1000 #Tf
        #ax0.plot(np.arange(0, T_final, 1), x_SF[0:T_final], color="red", label="Business as usual")
        ax0.plot(np.arange(0, T_final, 1), gmt[0:T_final], color="blue", label="Moderately incentivized clean energy")
        #ax0.plot(np.arange(0, T_final, 1), x_WF[0:T_final], color="green", label="Highly incentivized clean energy")

        ax0.set_xlabel("Time [yrs]")
        ax0.set_ylabel("$\Delta$ GMT [°C]")
        #ax0.set_ylim([0, 103])


        #sns.despine(bottom=True, left=True) #no right and upper border lines
        fig.tight_layout()
        fig.savefig("test_figures/scenarios_{}.png".format(nn))
        fig.savefig("test_figures/scenarios_{}.pdf".format(nn))
        # fig.show()
        fig.clf()
        plt.close()
        die
        """
        
        


        
        #necessary for break condition
        if len(output) != 0:
            #saving structure
            data = np.array(output, dtype=object)

            np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}.txt".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength), data[-1])
            time = data.T[0]
            state_gis = data.T[1]
            state_thc = data.T[2]
            state_wais = data.T[3]
            state_amaz = data.T[4]
            co2_series = data.T[-5] #co2 levels
            ecl_series = data.T[-4] #E_Cl
            ped_series = data.T[-2] #PED
            gmt_series = data.T[-1] #gmt series

            np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}_gmt.txt".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength), gmt_series)
            np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}_ped.txt".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength), ped_series)
            np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}_ecl.txt".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength), ecl_series)
            np.savetxt("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}_co2.txt".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength), co2_series)




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
            plt.savefig("{}/{}_feedbacks/network_{}_{}_{}/{}/feedbacks_care{:.2f}_{}_{:.2f}.pdf".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength))
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
            plt.savefig("{}/{}_feedbacks/network_{}_{}_{}/{}/gmtseries_care{:.2f}_{}_{:.2f}.pdf".format(long_save_name, namefile, 
                kk[0], kk[1], kk[2], str(mc_dir).zfill(4), care, scenario, strength))
            #plt.show()
            plt.clf()
            plt.close()
        


    
    # it is necessary to limit the amount of saved files
    # --> compose one pdf file for each network setting and remove the other time-files
    current_dir = os.getcwd()
    os.chdir("{}/{}_feedbacks/network_{}_{}_{}/{}/".format(long_save_name, namefile, kk[0], kk[1], kk[2], str(mc_dir).zfill(4)))
    pdfs = np.array(np.sort(glob.glob("feedbacks_care{:.2f}_{}_*.pdf".format(care, scenario)), axis=0))
    if len(pdfs) != 0.:
        merger = PdfFileMerger()
        for pdf in pdfs:
            merger.append(pdf)
        os.system("rm feedbacks_care{:.2f}_{}_*.pdf".format(care, scenario))
        merger.write("feedbacks_complete_care{:.2f}_{}.pdf".format(care, scenario))
        print("Complete PDFs merged - Part 1")

    pdfs = np.array(np.sort(glob.glob("gmtseries_care{:.2f}_{}_*.pdf".format(care, scenario)), axis=0))
    if len(pdfs) != 0.:
        merger = PdfFileMerger()
        for pdf in pdfs:
            merger.append(pdf)
        os.system("rm gmtseries_care{:.2f}_{}_*.pdf".format(care, scenario))
        merger.write("gmtseries_complete_care{:.2f}_{}.pdf".format(care, scenario))
        print("Complete PDFs merged - Part 2")
    os.chdir(current_dir)

print("Finish")
#end = time.time()
#print("Time elapsed until Finish: {}s".format(end - start))