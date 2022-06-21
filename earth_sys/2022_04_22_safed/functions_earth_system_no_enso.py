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
    def CUSPc(x1, x2, x):
        # only returns a value in case GMT is higher than lower boundary of the respective tipping element, otherwise return 0.0 (lower cap), N.B.: No upper cap
        if x >= x1:
            y1 = 0.0
            y2 = np.sqrt(4 / 27)
            y = (y2 - y1) / (x2 - x1) * (x - x1) + y1
            return y
        else:
            return 0.0


    def co2_gmt(co2):
        alpha = 4.841 # see IPCC, AR6, chapter 6, table A 6.2
        beta = 0.0906
        co2_preind = 280
        factor = 2.0226354355276586 #rescale the GMT such that 400ppm --> 1°C above pre-industrial levels
        gmt = (alpha*np.log(co2/co2_preind) + beta*(np.sqrt(co2)-np.sqrt(co2_preind)))/factor 
        return gmt


    #taking into account the uncertainty
    def unc_limits(ind_closest, switch_unc):
        if switch_unc == False:
            red_coop = 1.0
            return red_coop #returning 1.0 means that there is NO reduction in cooperativity with respect to uncertainty in the thresholds
        elif switch_unc == True:
            index = np.array([[0.8, 3.0], [1.4, 8.0], [1.0, 3.0], [2.0, 6.0]]) #limit of critical temperature thresholds of GIS, THC, WAIS, AMAZ
            unc = (index[ind_closest][1] - index[ind_closest][0])/(0.5*(index[ind_closest][0] + index[ind_closest][1]))

            #cooperation following Barrett, Dannenberg, 2014, Nature Climate Change (Fig.2, linear regression assumed if unc>=20%)
            P1 = np.array([2/15, 117.4/150.0]) 
            P2 = np.array([2/3, 77.2/150.0])

            #this is the reduced cooperativity factor taking into account the uncertainty in the thresholds
            red_coop = ((P2[1]-P1[1])/(P2[0]-P1[0])) * (unc-P1[0]) + P1[1]

            return red_coop
        else:
            print("No or a wrong switch_unc has been defined - please check!!")
            die