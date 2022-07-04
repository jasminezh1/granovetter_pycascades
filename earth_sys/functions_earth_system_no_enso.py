import numpy as np

"""
Here all global functions are stored - the functions are up to "choice": Here, linear functions are used
"""


class global_functions():
    """
    Linear feedback function to compute feedbacks. Maximal feedback is obtained from 2.0°C onwards.
    feedbacks for Arctic summer sea ice, mountain glaciers, GIS and WAIS are proven to be constant also for higher temperatures
    feedbacks for Amazon and feedbacks_steffen are computed for a temperature increase of 2.0°C until 2100, see paper
    """
    def feedback_function(GMT, fbmax):
        # fbmax = maximal feedback
        # only returns a value in case GMT is higher than lower boundary of the respective tipping element, otherwise return 0.0 (lower cap), N.B.: No upper cap
        if GMT <= 2.0:
            y = (fbmax / 2.0) * GMT
            return y
        elif GMT > 2.0:
            return fbmax
        else:
            raise Exception("GMT negativ: Feedbacks do not work for temperatures smaller 0!")

    # feedbacks of state dependent variables, linear increase of feedbacks between state -1 and +1 for GIS, WAIS, THC, NINO and AMAZ
    def state_feedback(state, fbmax):
        if state >= -1 and state <= 1:
            y = fbmax / 2 * (state + 1)
        elif state < -1:
            y = 0.
        elif state > +1:
            y = fbmax
        return y

    # c = c(GMT) where tipping occurs at sqrt(4/27) ~ 0.38
    # Linear function through two points maps GMT --> c, where x-values represent GMT and y-values represent CUSP-c values
    # function rescales the critical gmt? tipping element tips at sqrt 4/27
    def CUSPc(x1, x2, x):
        # only returns a value in case GMT is higher than lower boundary of the respective tipping element, otherwise return 0.0 (lower cap), N.B.: No upper cap
        if x >= x1:
            y1 = 0.0
            y2 = np.sqrt(4 / 27)
            y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
            #print("VALUES ", x1, " ", x2, " ", x, " ", y)
            return y
        else:
            return 0.0


    def co2_gmt(co2):
        alpha = 4.841 # see IPCC, AR6, chapter 6, table A 6.2, second formula
        beta = 0.0906
        co2_preind = 280
        factor = 2.0226354355276586 #rescale the GMT such that 400ppm --> 1°C above pre-industrial levels
        gmt = (alpha*np.log(co2/co2_preind) + beta*(np.sqrt(co2)-np.sqrt(co2_preind)))/factor 
        return gmt

    def gmt_co2(gmt): # see IPCC, AR6, chapter 6, table A 6.2, first formula solved for co2
        alpha = 5.35
        co2_preind = 280
        factor = 0.8438681 #rescale the GMT such that 400ppm --> 1°C above pre-industrial levels
        co2 = co2_preind*np.exp(gmt/alpha)/factor
        return co2


    #taking into account the uncertainty based on Barrett, Dannenber papers
    def unc_limits_barrett(ind_closest, unc_temp):
        unc = (unc_temp[ind_closest][1] - unc_temp[ind_closest][0])/(0.5*(unc_temp[ind_closest][0] + unc_temp[ind_closest][1]))

        #cooperation following Barrett, Dannenberg, 2014, Nature Climate Change (Fig.2, linear regression assumed if unc>=20%)
        P1 = np.array([2/15, 117.4/150.0]) 
        P2 = np.array([2/3, 77.2/150.0])

        #this is the reduced cooperativity factor taking into account the uncertainty in the thresholds
        red_coop = ((P2[1]-P1[1])/(P2[0]-P1[0])) * (unc-P1[0]) + P1[1]
        return red_coop


    #defined as the uncertainty of thresholds
    def unc_limits(ind_closest, unc_temp):
        kappa_know_unc = np.random.uniform(unc_temp[ind_closest][0], unc_temp[ind_closest][1], 1) #randomly chose between the limits
        return kappa_know_unc
 