from __future__ import print_function
import numpy as np
import sys
import math
import mpmath as mpm
from datetime import datetime
import multiprocessing as mp
import warnings

from backscatter import dmap as dm
from backscatter import config
from backscatter import hdw_info as hdw


MIN_LAGS = int(config.get("fitacf","minimum_lags"))
ACF_SNR_CUTOFF = float(config.get("fitacf","acf_snr_cutoff"))
ALPHA_CUTOFF = float(config.get("fitacf","alpha_cutoff"))
FLUCTUATION_CUTOFF_COEFF = int(config.get("fitacf","fluctuation_cutoff_coeff"))
FITACF_REVISION_MAJOR = int(config.get("fitacf","fitacf_revision_major"))
FITACF_REVISON_MINOR = int(config.get("fitacf","fitacf_revision_minor"))
V_MAX = float(config.get("fitacf","v_max"))
W_MAX = float(config.get("fitacf","w_max"))


C = 299792458.0


#for printing to stderr easily and portably
def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def create_lag_list(raw_data):
    """Creates a list of lag dictionaries from raw data
    
    This method uses the mplgs, ptab, mppul, ltab, mpinc,
    and smsep fields of the raw data to create a dictionary
    for each lag. Each lag dictionary contains a field for
    it's number, the pulses used to make the lag, the indices
    at which those pulses are located in ptab, and the sample
    bases.
    
    Lag fields {'number','pulses','pulse2_idx','pulse1_idx','sample_base1','sample_base2'}
    :param raw_data: a dictionary of raw data parameters
    :returns: list of lag dictionaries

    """

    mplgs = raw_data['mplgs']
    pulses = raw_data['ptab']
    mppul = raw_data['mppul']
    lag_table = raw_data['ltab']
    mpinc = raw_data['mpinc']
    smsep = raw_data['smsep']

    lags = []
    for i in range(0,mplgs):
        lag = {}
        lag['number'] = lag_table[i][1] - lag_table[i][0]

        for idx,pulse in enumerate(pulses):
            if lag_table[i][1] == pulse:
                lag['pulse2_idx'] = idx

            if lag_table[i][0] == pulse:
                lag['pulse1_idx'] = idx 

        lag['sample_base1'] = lag_table[i][0] * (mpinc/smsep)
        lag['sample_base2'] = lag_table[i][1] * (mpinc/smsep)
        lag['pulses'] = lag_table[i]
        lags.append(lag)

    return lags

class PowerDataPoints(object):
    """Contains the power data points for a particular range

    This class contains creates an array of log powers, sigmas,
    and t values for a range created from the raw data. 

    """

    def __init__(self, raw_data,lags,range_obj):
        self.log_pwrs = None
        self.sigmas = None
        self.t = None

        self.create_arrays(raw_data,lags,range_obj)


    def create_arrays(self,raw_data,lags,range_obj):
        """Creates the data arrays associated with a range

        From the raw data, the magnitude of the power is found.
        It is then normalized for the calculation of sigma. 
        After sigma is found, t is determined by multiplying 
        lag numbers by the fundamental spacing.

        :param raw_data: a dictionary of raw data parameters
        :param lags: list of lag dictionaries
        :param range_obj: The range object this data is to associated with

        """
        acfd = raw_data['acfd'][range_obj.range_number]
        mplgs = raw_data['mplgs']
        pwr0 = raw_data['pwr0'][range_obj.range_number]
        nave = raw_data['nave']
        mpinc = raw_data['mpinc']

        real = acfd[:,0]
        imag = acfd[:,1]
        real_2 = real**2
        imag_2 = imag**2

        pwrs = np.sqrt(np.add(real_2,imag_2))

        pwr_normalized = np.divide(pwrs,pwr0)
        pwr_normalized_2 = pwr_normalized**2
        inverse_alpha_2 = np.reciprocal(range_obj.alpha_2)

        sigmas = [pwr0 * np.sqrt((pn_2 + ia_2)/(2 * nave)) for (pwr, pn_2, ia_2) in zip(pwrs, pwr_normalized_2,inverse_alpha_2)]

        self.sigmas = np.array(sigmas)

        t_values = [lag['number'] * mpinc * 1.0e-6 for (pwr, lag) in zip(pwrs,lags)]

        self.t = np.array(t_values)
        
        #there will for sure be log of 0 here, but we filter it and 
        #dont need to be warned
        warnings.simplefilter("ignore")
        self.log_pwrs = np.log(pwrs)




    def remove_bad_points(self,bad_indices):
        """Removes data points that are to be excluded from fitting

        :param bad_indices: a list of indices of points to remove

        """

        if not bad_indices:
            return

        num_points = len(self.log_pwrs)

        #Removes any indices that are larger than length of data array
        #incase those points were removed already

        if max(bad_indices) >= num_points:
            bad_indices = [bi for bi in bad_indices if bi < num_points]

        mask = np.ones(num_points, np.bool)
        mask[bad_indices] = 0

        self.log_pwrs = self.log_pwrs[mask]
        self.sigmas = self.sigmas[mask]
        self.t = self.t[mask]
        # self.log_pwrs = np.delete(self.log_pwrs,bad_indices)
        # self.sigmas = np.delete(self.sigmas,bad_indices)
        # self.t = np.delete(self.t,bad_indices)

    def remove_inf_points(self,non_inf_indices):

        self.log_pwrs = self.log_pwrs[non_inf_indices]
        self.sigmas = self.sigmas[non_inf_indices]
        self.t = self.t[non_inf_indices]


class PhaseDataPoints(object):
    """Contains phase data points for a particular range. 

    Phase data points can apply to both ACF or XCF phase. This 
    class is used for both velocity and elevation points. This 
    class creates an array of phases, placeholder sigmas(alpha_2),
    and t values for a range created from the raw data. 

    """
    
    def __init__(self,raw_data,phase_type,lags,range_obj):
        self.phases = None
        self.sigmas = None
        self.t = None

        self.create_arrays(raw_data,phase_type,lags,range_obj)


    def create_arrays(self,raw_data,phase_type,lags,range_obj):
        """Creates the data arrays associated with a range

        From the raw data, phase is determined for ACF or XCF data.
        Sigmas are determined after power fitting, so alpha_2 is
        used a placeholder at this point. After sigma is found, t 
        is determined by multiplying lag numbers by the fundamental 
        spacing.

        :param raw_data: a dictionary of raw data parameters
        :param phase_type: "acfd" or "xcfd" to select which data arrays to use
        :param lags: list of lag dictionaries
        :param range_obj: The range object this data is to associated with

        """        

        acfd = raw_data[phase_type][range_obj.range_number]
        mplgs = raw_data['mplgs']
        mpinc = raw_data['mpinc']

        real = acfd[:,0]
        imag = acfd[:,1]

        self.phases = np.arctan2(imag,real)

        self.sigmas = np.copy(range_obj.alpha_2)

        self.t = np.array([lag['number'] * mpinc * 1.0e-6 for lag in lags])

    def remove_bad_points(self,bad_indices):
        """Removes data points that are to be excluded from fitting

        :param bad_indices: a list of indices of points to remove

        """
        if not bad_indices:
            return
          
        num_points = len(self.phases)

        #Removes any indices that are larger than length of data array
        #incase those points were removed already
        if max(bad_indices) >= num_points:
            bad_indices = [bi for bi in bad_indices if bi < num_points]

        mask = np.ones(num_points, np.bool)
        mask[bad_indices] = 0
        self.phases = self.phases[mask]
        self.sigmas = self.sigmas[mask]
        self.t = self.t[mask]
        # self.phases = np.delete(self.phases,bad_indices)
        # self.sigmas = np.delete(self.sigmas,bad_indices)
        # self.t = np.delete(self.t,bad_indices)

    def set_sigmas(self,sigmas):
        """Reassign sigma values

        :param sigmas: an array of new sigmas

        """
        self.sigmas = sigmas

    def set_phases(self,phases):
        """Reassign phases

        :param phases: an array of new phases

        """
        self.phases = phases


class Range(object):
    """This class holds all the data associated with a range to be fit

    The Range class extracts what is necessary from the raw data for
    a particular range in order to prepare for a fit. This class computes
    the cross-range interference for a range and then generates the alpha_2
    values for each lag. Phases, elevations, and power data points are then
    constructed and calculated from the raw data.


    """

    def __init__(self, range_number,raw_data,lags):
        self.range_number = range_number
        self.CRI = self.find_cri(raw_data)
        self.alpha_2 = self.find_alphas(raw_data,lags)
        self.phases = PhaseDataPoints(raw_data,'acfd',lags,self)
        self.elevs = PhaseDataPoints(raw_data,'xcfd',lags,self)
        self.pwrs = PowerDataPoints(raw_data,lags,self)

        self.linear_pwr_fit = None
        self.quadratic_pwr_fit = None
        self.linear_pwr_fit_err = None
        self.quadratic_pwr_fit_err = None
        self.phase_fit = None
        self.elev_fit = None

    def find_cri(self,raw_data):
        """Creates an array of cross range interference for each pulse
        
        :param raw_data: a dictionary of raw data parameters

        """

        smsep = raw_data['smsep']
        mpinc = raw_data['mpinc']
        txpl = raw_data['txpl']
        mppul = raw_data['mppul']
        pulses = raw_data['ptab']
        nrang = raw_data['nrang']
        pwr0 = raw_data['pwr0']

        if smsep != 0:
            tau = mpinc/smsep
        else:
            eprint("r_overlap: WARNING, using txpl instead of smsep...\n")
            tau = mpinc/txpl

        cri_for_pulses = np.ones(mppul)
        for pulse_to_check in range(0,mppul):
            total_cri = 0.0

            for pulse in range(0,mppul):
                pulse_diff = pulses[pulse_to_check] - pulses[pulse]
                range_to_check = pulse_diff * tau + self.range_number

                if (pulse != pulse_to_check and
                    0 <= range_to_check and
                    range_to_check < nrang):

                    total_cri = total_cri + pwr0[range_to_check]

            cri_for_pulses[pulse_to_check] = total_cri

        return cri_for_pulses

    def find_alphas(self,raw_data,lags):
        """From cross-range interference, computes alpha_2 for each lag

        :param raw_data: a dictionary of raw data parameters
        :param lags: a list of lag dictionaries
        :returns: an array of alphas for each lag

        """

        pwr0 = raw_data['pwr0']

        alpha_2 = np.ones(len(lags))
        for idx,lag in enumerate(lags):
            pulse1_cri = self.CRI[lag['pulse1_idx']]
            pulse2_cri = self.CRI[lag['pulse2_idx']]

            lag_0_pwr = pwr0[self.range_number]
            alpha_2[idx] = lag_0_pwr**2/((lag_0_pwr + pulse1_cri) * (lag_0_pwr + pulse2_cri))

        return alpha_2

    def remove_bad_alphas(self,bad_indices):
        """Remove alpha_2 that are associated with bad data points

        :param bad_indices: a list of indices of points to remove

        """
        if not bad_indices:
            return
          
        num_points = len(self.alpha_2)

        #Removes any indices that are larger than length of data array
        #incase those points were removed already
        if max(bad_indices) >= num_points:
            bad_indices = [bi for bi in bad_indices if bi < num_points]


        mask = np.ones(len(self.alpha_2), np.bool)
        mask[bad_indices] = 0
        self.alpha_2 = self.alpha_2[mask]


                                                
class NoisePower(object):
    """This class holds static methods used to calculate noise

    """

    @staticmethod   
    def cutoff_power_correction(raw_data):
        """Computes the correction factor for noise power

        Computes a correction factor for the noise. Without this
        factor, noise is underestimated by 30%-40%.

        :param raw_data: a dictionary of raw data parameters
        :returns: power correction factor

        """
        nave = raw_data['nave']
        nrang = raw_data['nrang']
    
        std_dev = 1.0/math.sqrt(nave)
    
        i = 0
        cumulative_pdf = 0.0
        cumulative_pdf_x_norm_pwr = 0
        while cumulative_pdf < (10.0/nrang):
            #Normalised power for calculating model PDF (Gaussian)
            normalized_pwr = i/1000.0
    
            x  = -((normalized_pwr - 1.0)**2/(2.0 * std_dev**2))
            pdf = math.exp(x)/std_dev/math.sqrt(2 * math.pi)/1000
            cumulative_pdf = cumulative_pdf + pdf
    
            #Cumulative value of PDF*x -- needed for calculating the mean
            cumulative_pdf_x_norm_pwr = cumulative_pdf_x_norm_pwr + pdf * normalized_pwr
    
            i = i + 1
    
        #Correcting factor as the inverse of a normalised mean
        corr = 1.0/(cumulative_pdf_x_norm_pwr/cumulative_pdf)
        return corr

    @staticmethod
    def acf_cutoff_pwr(raw_data):
        """Determines the flucuation level for which ACFs are pure noise

        Uses the ten weakest ACFS to determine noise level. A noise
        correction is applied to reduce underestimation.

        :param raw_data: a dictionary of raw data parameters
        :returns: estimate of ACF fluctuation level 

        """
        sorted_pwr_levels = np.sort(raw_data['pwr0'])
        
        i, j = 0, 0
        min_pwr = 0
        nrang = raw_data['nrang']
        while (j < 10 and i < nrang/3):
            if sorted_pwr_levels[i] > 0.0:
                j = j + 1
            min_pwr = min_pwr + sorted_pwr_levels[i]
            i = i + 1
    
        if j <= 0:
            j = 1
    
        min_pwr = min_pwr/j * NoisePower.cutoff_power_correction(raw_data)
    
        if min_pwr < 1.0:
            min_pwr = raw_data['noise.search']
    
        return min_pwr

class LeastSquaresValues(object):
    """This class simply holds all the least squares values associated
    with the fitting algorithm outlined in NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING

    """
    def __init__(self):
        self.S = 0.0
        self.S_x = 0.0
        self.S_y = 0.0
        self.S_xx = 0.0
        self.S_xy = 0.0
        self.delta = 0.0
        self.a = 0.0
        self.b = 0.0
        self.sigma_2_a = 0.0
        self.sigma_2_b = 0.0
        self.delta_a = 0.0
        self.delta_b = 0.0
        self.cov_ab = 0.0
        self.r_ab = 0.0
        self.Q = 0.0
        self.chi_2 = 0.0

    
class LeastSquaresFitting(object):
    """This class holds all the methods needed to fit rawacf data

    This class holds methods for 1 parameter straight line fits,
    2 parameter straight line fits, and 2 parameter quadratic fits.

    """

    def __init__(self,confidence,DoF):
        self.delta_chi_2 = [[1.00,2.30],
                            [2.71,4.61],
                            [4.00,6.17],
                            [6.63,9.21],
                            [9.00,11.8],
                            [15.1,18.4]]
        self.confidence = confidence - 1
        self.DoF = DoF - 1



    def find_sums(self,lst_sqrs,x_arr,y_arr,sigmas,fit_type):
        """Computes the sums needed for linear least squares equations

        :param lst_sqrs: LeastSquaresValues object to fill
        :param x_arr: array with x-axis data
        :param y_arr: array with y-axis data
        :param sigmas: weighting for y-axis data
        :param fit_type: selects between 'linear' and 'quadratic' fits

        """
        sigma_nonzero = np.nonzero(sigmas)

        x_arr = x_arr[sigma_nonzero]
        y_arr = y_arr[sigma_nonzero]
        sigmas = sigmas[sigma_nonzero]

        sigma_2 = sigmas * sigmas

        if fit_type == 'linear':
            lst_sqrs.S = np.sum(np.reciprocal(sigma_2))
            lst_sqrs.S_x = np.sum(x_arr/sigma_2)
            lst_sqrs.S_y = np.sum(y_arr/sigma_2)
            lst_sqrs.S_xx = np.sum((x_arr * x_arr)/sigma_2)
            lst_sqrs.S_xy = np.sum((x_arr * y_arr)/sigma_2) 
        elif fit_type == 'quadratic':
            x_2 = x_arr * x_arr
            lst_sqrs.S = np.sum(np.reciprocal(sigma_2)) 
            lst_sqrs.S_x = np.sum(x_2/sigma_2)
            lst_sqrs.S_y = np.sum(y_arr/sigma_2)
            lst_sqrs.S_xx = np.sum((x_2 * x_2)/sigma_2)
            lst_sqrs.S_xy = np.sum((x_2 * y_arr)/sigma_2)
        else:
            error_msg = "Invalid fit type {0} in find_sums()".format(fit_type)
            raise ValueError(error_msg)
            

    def one_parameter_line_fit(self,x_arr,y_arr,sigmas,num_points):
        """Computes a fit for the model y = bx

        :param x_arr: array with x-axis data
        :param y_arr: array with y-axis data
        :param sigmas: weighting for y-axis data
        :param num_points: number of data points
        :returns: LeastSquaresValues with computed values

        """
        lst_sqrs = LeastSquaresValues()

        self.find_sums(lst_sqrs,x_arr,y_arr,sigmas,'linear')

        S_xx = lst_sqrs.S_xx
        S_xy = lst_sqrs.S_xy

        lst_sqrs.delta = 0.0
        lst_sqrs.a = 0.0
        lst_sqrs.b = S_xy / S_xx
        lst_sqrs.sigma_2_a = 0.0
        lst_sqrs.sigma_2_b = 1/S_xx
        lst_sqrs.cov_ab = 0.0
        lst_sqrs.r_ab = 0.0

        dc2 = self.delta_chi_2[self.confidence][self.DoF]
        lst_sqrs.delta_a = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_a)
        lst_sqrs.delta_b = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_b)


        self.find_chi_2(lst_sqrs,x_arr,y_arr,sigmas,'linear')

        #lst_sqrs.Q = mpm.gammainc((num_points - 2) * .5, lst_sqrs.chi_2 * 0.5)

        return lst_sqrs

    def two_parameter_line_fit(self,x_arr,y_arr,sigmas,num_points):
        """Computes a fit for the model y = bx + a

        :param x_arr: array with x-axis data
        :param y_arr: array with y-axis data
        :param sigmas: weighting for y-axis data
        :param num_points: number of data points
        :returns: LeastSquaresValues with computed values

        """        
        lst_sqrs = LeastSquaresValues()

        self.find_sums(lst_sqrs,x_arr,y_arr,sigmas,'linear')

        S = lst_sqrs.S
        S_x = lst_sqrs.S_x
        S_y = lst_sqrs.S_y
        S_xx = lst_sqrs.S_xx
        S_xy = lst_sqrs.S_xy

        lst_sqrs.delta = S * S_xx - S_x * S_x
        lst_sqrs.a = (S_xx * S_y - S_x * S_xy)/lst_sqrs.delta
        lst_sqrs.b = (S * S_xy - S_x * S_y)/lst_sqrs.delta
        lst_sqrs.sigma_2_a = S_xx/lst_sqrs.delta
        lst_sqrs.sigma_2_b = S/lst_sqrs.delta
        lst_sqrs.cov_ab = -1 * S_x/lst_sqrs.delta
        lst_sqrs.r_ab = -1 * S_x/math.sqrt(S * S_xx)

        dc2 = self.delta_chi_2[self.confidence][self.DoF]
        lst_sqrs.delta_a = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_a)
        lst_sqrs.delta_b = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_b)

        self.find_chi_2(lst_sqrs,x_arr,y_arr,sigmas,'linear')

        #lst_sqrs.Q = mpm.gammainc((num_points - 2) * .5, lst_sqrs.chi_2 * 0.5)

        return lst_sqrs

    def quadratic_fit(self,x_arr,y_arr,sigmas,num_points):
        """Computes a fit for the model y = bx^2 + a

        :param x_arr: array with x-axis data
        :param y_arr: array with y-axis data
        :param sigmas: weighting for y-axis data
        :param num_points: number of data points
        :returns: LeastSquaresValues with computed values

        """ 
        lst_sqrs = LeastSquaresValues()

        self.find_sums(lst_sqrs,x_arr,y_arr,sigmas,'quadratic')

        S = lst_sqrs.S
        S_x = lst_sqrs.S_x
        S_y = lst_sqrs.S_y
        S_xx = lst_sqrs.S_xx
        S_xy = lst_sqrs.S_xy

        lst_sqrs.delta = S * S_xx - S_x**2
        lst_sqrs.a = (S_xx * S_y - S_x * S_xy)/lst_sqrs.delta
        lst_sqrs.b = (S * S_xy - S_x * S_y)/lst_sqrs.delta
        lst_sqrs.sigma_2_a = S_xx/lst_sqrs.delta
        lst_sqrs.sigma_2_b = S/lst_sqrs.delta
        lst_sqrs.cov_ab = -1 * S_x/lst_sqrs.delta
        lst_sqrs.r_ab = -1 * S_x/math.sqrt(S * S_xx)

        dc2 = self.delta_chi_2[self.confidence][self.DoF]
        lst_sqrs.delta_a = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_a)
        lst_sqrs.delta_b = math.sqrt(dc2) * math.sqrt(lst_sqrs.sigma_2_b)

        self.find_chi_2(lst_sqrs,x_arr,y_arr,sigmas,'quadratic')

        #lst_sqrs.Q = mpm.gammainc((num_points - 2) * .5, lst_sqrs.chi_2 * 0.5)

        return lst_sqrs


    def find_chi_2(self,lst_sqrs,x_arr,y_arr,sigmas,fit_type):
        """Computes the chi-square statistic of the fit

        :param lst_sqrs: LeastSquaresValues object with parameters a and b
        :param x_arr: array with x-axis data
        :param y_arr: array with y-axis data
        :param sigmas: weighting for y-axis data
        :param fit_type: selects between 'linear' and 'quadratic' fits

        """ 
        sigma_nonzero = np.nonzero(sigmas)

        x_arr = x_arr[sigma_nonzero]
        y_arr = y_arr[sigma_nonzero]
        sigmas = sigmas[sigma_nonzero]

        if fit_type == 'linear':
            chi = ((y_arr - lst_sqrs.a) - (lst_sqrs.b * x_arr)) / sigmas
            lst_sqrs.chi_2 = np.sum(chi**2)
        elif fit_type == 'quadratic':
            chi = ((y_arr - lst_sqrs.a) - (lst_sqrs.b * (x_arr**2))) / sigmas
            lst_sqrs.chi_2 = np.sum(chi**2)


class ACFFitting(object):
    """This class is a container for methods that initiate fitting
    to the data

    There are static methods to perform fits for power, velocity and
    elevation, as well the methods to calculate sigmas and unwrap phase
    for ACF and XCF phases.
    """

    @staticmethod
    def acf_pwr_fitting(range_list):
        """Initiates least square fitting for power

        Performs a least squares fit for two parameter linear and              
        quadratic models to ACF power for each range. Performs additional
        linear and quadratic fits using log-corrected sigmas for correct 
        errors. 

        :param range_list: list of Range objects with data to fit

        """
        lst_sqrs_fit = LeastSquaresFitting(1,1)
        two_param_line_fit = lst_sqrs_fit.two_parameter_line_fit
        quad_fit = lst_sqrs_fit.quadratic_fit

        for range_obj in range_list:
            log_pwrs = range_obj.pwrs.log_pwrs
            sigmas = range_obj.pwrs.sigmas
            t = range_obj.pwrs.t

            if len(log_pwrs) != len(sigmas) != len(t):
                error_msg = """The length of the data point arrays dont agree in 
                power fitting! {0} log_pwr points, {1} sigma points, {2} t points.
                """.format(len(log_pwrs),len(sigmas),len(t))
                raise ValueError(error_msg)
                
            else:
                num_points = len(log_pwrs)

            range_obj.linear_pwr_fit = two_param_line_fit(t,log_pwrs,sigmas,num_points)
            range_obj.quadratic_pwr_fit = quad_fit(t,log_pwrs,sigmas,num_points)

            log_corrected_sigmas = sigmas / np.exp(log_pwrs)

            range_obj.linear_pwr_fit_err = two_param_line_fit(t,log_pwrs,log_corrected_sigmas,num_points)
            range_obj.quadratic_pwr_fit_err = quad_fit(t,log_pwrs,log_corrected_sigmas,num_points)

    @staticmethod
    def acf_phase_fitting(range_list):
        """Initiates least square fitting for ACF phase

        Performs a least squares fit for the one parameter linear model to 
        ACF phase for each range.

        :param range_list: list of Range objects with data to fit

        """
        lst_sqrs_fit = LeastSquaresFitting(1,1)
        one_param_line_fit = lst_sqrs_fit.one_parameter_line_fit

        for range_obj in range_list:
            phase_values = range_obj.phases.phases
            sigma_values = range_obj.phases.sigmas
            t_values = range_obj.phases.t

            if len(phase_values) != len(sigma_values) != len(t_values):
                error_msg = """The length of the data point arrays dont agree in 
                phase fitting! {0} phase points, {1} sigma points, {2} t points.
                """.format(len(phase_values),len(sigmas),len(t))
                raise ValueError(error_msg)
            else:
                num_points = len(phase_values)

            range_obj.phase_fit = one_param_line_fit(t_values,phase_values,sigma_values,num_points)

    @staticmethod
    def xcf_phase_fitting(range_list):
        """Initiates least square fitting for XCF phase

        Performs a least squares fit for the two parameter linear model to 
        XCF phase for each range.

        :param range_list: list of Range objects with data to fit

        """
        lst_sqrs_fit = LeastSquaresFitting(1,1)
        two_param_line_fit = lst_sqrs_fit.two_parameter_line_fit

        for range_obj in range_list:
            elev_values = range_obj.elevs.phases
            sigma_values = range_obj.elevs.sigmas
            t_values = range_obj.elevs.t

            if len(elev_values) != len(sigma_values) != len(t_values):
                error_msg = """The length of the data point arrays dont agree in 
                elevation fitting! {0} phase points, {1} sigma points, {2} t points.
                """.format(len(elev_values),len(sigmas),len(t))
                raise ValueError(error_msg)
            else:
                num_points = len(elev_values)
                
            range_obj.elev_fit = two_param_line_fit(t_values,elev_values,sigma_values,num_points)


    @staticmethod
    def calculate_phase_and_elev_sigmas(range_list,raw_data):
        """Calculates correct weightings for ACF and XCF phases

        ACF phase and XCF phase sigmas can only be computed after
        ACF power has been fitted. This method computes the correct
        values for sigmas for each range.

        :param range_list: list of Range objects with phase data and fitted power
        :param raw_data: a dictionary of raw data parameters

        """
        nave = raw_data['nave']

        for range_obj in range_list:
            phases = range_obj.phases
            pwrs = range_obj.pwrs
            elevs = range_obj.elevs

            phase_inverse_alpha_2 = 1/phases.sigmas
            elev_inverse_alpha_2 = 1/elevs.sigmas

            #phase and elevation have same t values
            pwr_values = np.exp(-1 * math.fabs(range_obj.linear_pwr_fit.b) * phases.t)
            inverse_pwr_2_values = 1/(pwr_values**2)
            
            #print(len(elevs.sigmas),len(elev_inverse_alpha_2),len(phases.sigmas),len(pwr_values))

            phase_numerator = ((phase_inverse_alpha_2 * inverse_pwr_2_values) - 1)
            elev_numerator = ((elev_inverse_alpha_2 * inverse_pwr_2_values) - 1)
            denominator = 2 * nave

            phase_sigmas = np.sqrt((phase_numerator/denominator))
            elev_sigmas = np.sqrt((elev_numerator/denominator))
            
            if np.isnan(phase_sigmas).any() or np.isinf(phase_sigmas).any():
                error_string = "Phase sigmas bad at range {0} -- phase_inverse_alphas,pwr_values"
                error_string.format(range_obj.range_number)
                eprint(error_string,phase_inverse_alpha_2,pwr_values)

            if np.isnan(elev_sigmas).any() or np.isinf(elev_sigmas).any():
                error_string = "Elevation sigmas bad at range {0} -- elev_inverse_alphas,pwr_values"
                error_string.format(range_obj.range_number)
                eprint(error_string,elev_inverse_alpha_2,pwr_values)

            phase_sigmas = np.array([math.pi if sigma > math.pi else sigma for sigma in phase_sigmas])
            elev_sigmas = np.array([math.pi if sigma > math.pi else sigma for sigma in elev_sigmas])

            """Since lag 0 phase is included for elevation fit, we set lag 0 sigma the
            same as lag 1 sigma"""
            elev_sigmas[0] = elev_sigmas[1]

            phases.set_sigmas(phase_sigmas)
            elevs.set_sigmas(elev_sigmas)

    @staticmethod
    def acf_phase_unwrap(range_list):
        """Unwraps the ACF phase data to be able to fit a straight line

        Takes phase data in a wrapping domain and converts it to a sloped
        line using a 2 stage iterative process. This is to be able to fit 
        using linear least squares. Performed after sigma values are calculated.

        :param range_list: list of Range object with phase data

        """
        phase_correction = ACFFitting.phase_correction

        for range_obj in range_list:
            phases = range_obj.phases

            phase_values = phases.phases
            sigma_values = phases.sigmas
            t_values = phases.t

            phase_prev = phase_values[0]
            sigma_prev = sigma_values[0]
            t_prev = t_values[0]

            
            slope_numerator, slope_denominator = 0, 0

            #This is to skip the first element in lists
            iterator = np.nditer([phase_values,sigma_values,t_values])
            iterator.iternext()
            while not iterator.finished:
                phase_curr = iterator[0]
                sigma_curr = iterator[1]
                t_curr = iterator[2]

                phase_diff = phase_curr - phase_prev

                sigma_bar = (sigma_curr + sigma_prev)/2

                t_diff = t_curr - t_prev

                if np.fabs(phase_diff) < math.pi:
                    slope_numerator = slope_numerator + phase_diff/(sigma_bar**2)/t_diff
                    slope_denominator = slope_denominator + 1/(sigma_bar**2)

                phase_prev = phase_curr
                sigma_prev = sigma_curr
                t_prev = t_curr

                iterator.iternext()


            slope_estimate = slope_numerator / slope_denominator
            new_phases = np.array([phase_correction(slope_estimate,phase,t) for phase,t in zip(phase_values,t_values)])

            iterator = np.nditer([new_phases,sigma_values,t_values])

            S_xx, S_xy = 0.0, 0.0
            while not iterator.finished:
                phase = iterator[0]
                sigma = iterator[1]
                t = iterator[2]
                if sigma > 0.0:
                    S_xy = S_xy + (phase * t)/(sigma**2)
                    S_xx = S_xx + (t**2)/(sigma**2)

                iterator.iternext()

            slope_estimate = S_xy / S_xx
            new_phases = np.array([phase_correction(slope_estimate,phase,t) for phase,t in zip(new_phases,t_values)])

            range_obj.phases.set_phases(new_phases)
    
    @staticmethod
    def xcf_phase_unwrap(range_list):
        """Unwraps the ACF phase data to be able to fit a straight line

        Takes phase data in a wrapping domain and converts it to a sloped
        line. XCF unwrapping only needs 1 stage of iteration because it uses 
        the ACF phase fit as an initial guess. This is to be able to fit 
        using linear least squares. Performed after ACF phase is fit.

        :param range_list: list of Range object with phase data

        """
        phase_correction = ACFFitting.phase_correction

        for range_obj in range_list:
            elevs = range_obj.elevs

            elev_phase_values = elevs.phases
            sigma_values = elevs.sigmas
            t_values = elevs.t

            if range_obj.phase_fit is None:
                error_msg = """Phase fit must be defined in order to begin
                XCF phase unwrap!"""
                raise ValueError(error_msg)

            new_phases = np.array([phase_correction(range_obj.phase_fit.b,elev_phase,t) for elev_phase,t in zip(elev_phase_values,t_values)])

            iterator = np.nditer([new_phases,sigma_values,t_values])
            S_xx, S_xy = 0.0, 0.0
            while not iterator.finished:
                elev_phase = iterator[0]
                sigma = iterator[1]
                t = iterator[2]
                if sigma > 0.0:
                    S_xy = S_xy + (elev_phase * t)/(sigma**2)
                    S_xx = S_xx + (t**2)/(sigma**2)

                iterator.iternext()

            slope_estimate = S_xy / S_xx
            new_phases = np.array([phase_correction(slope_estimate,elev_phase,t) for elev_phase,t in zip(new_phases,t_values)])

            range_obj.elevs.set_phases(new_phases)

    @staticmethod
    def phase_correction(slope_estimate,phase_value,t_value):
        """Adds the estimated number of 2*pi corrections to phase values

        :param slope_estimate: predicted slope from iterative unwrap
        :param phase_value: a particular phase to correct
        :param t_value: a value in time for a particular phase point
        :returns: phase shifted by number 2*pi corrections

        """

        phase_predicted = slope_estimate * t_value
        phase_correction = round((phase_predicted - phase_value)/(2 * math.pi))

        return phase_value + phase_correction * 2 * math.pi



class Filtering(object):
    """This class contains static methods relating to the filtering
    of bad data from the consideration of fitting.

    """

    @staticmethod
    def mark_bad_samples(raw_data):
        """Mark the samples that are blacked out by TX overlap

        :param raw_data: a dictionary of raw data parameters

        """
        mppul = raw_data['mppul']
        pulses = raw_data['ptab']
        mpinc = raw_data['mpinc']
        txpl = raw_data['txpl']
        smsep = raw_data['smsep']
        lagfr = raw_data['lagfr']

        i = -1
        ts, t1, t2, sample = lagfr, 0, 0, 0
        bad_samples = []
        while i < (mppul -1):
            # first, skip over any pulses that occur before the first sample
            # defines transmitter pulse blanket window
            while ((ts > t2) and (i < (mppul - 1))):
                i = i + 1
                t1 = pulses[i] * mpinc - txpl/2
                t2 = t1 + 3 * txpl/2 + 100

            # we now have a pulse that occurs after the current sample.  Start
            # incrementing the sample number until we find a sample that lies
            # within the pulse
            while (ts < t1):
                sample = sample + 1
                ts = ts + smsep

            # ok, we now have a sample which occurs after the pulse starts.
            # check to see if it occurs before the pulse ends, and if so, mark
            # it as a bad sample
            while ((ts >= t1) and (ts <= t2)):
                bad_samples.append(sample)
                sample = sample + 1
                ts = ts + smsep

        return bad_samples

    @staticmethod
    def filter_tx_overlapped_lags(raw_data,lags,range_list):
        """Remove data points affected by TX overlapped lags

        :param raw_data: a dictionary of raw data parameters
        :param lags: list of lag dictionaries
        :param range_list: A list of Range objects with data points

        """

        bad_samples = Filtering.mark_bad_samples(raw_data)

        for range_obj in range_list:

            bad_indices = []
            for idx,lag in enumerate(lags):
                sample1 = lag['sample_base1'] + range_obj.range_number
                sample2 = lag['sample_base2'] + range_obj.range_number

                if (sample1 in bad_samples) or (sample2 in bad_samples):
                    bad_indices.append(idx)

            range_obj.pwrs.remove_bad_points(bad_indices)
            range_obj.phases.remove_bad_points(bad_indices)
            range_obj.elevs.remove_bad_points(bad_indices)
            range_obj.remove_bad_alphas(bad_indices)



    @staticmethod
    def filter_inf_lags(range_list):

        for range_obj in range_list:
            log_pwrs = range_obj.pwrs.log_pwrs

            inf_indices = [idx for idx,log_pwr in enumerate(log_pwrs) if not np.isfinite(log_pwr)]
 
            range_obj.pwrs.remove_bad_points(inf_indices)
            range_obj.remove_bad_alphas(inf_indices)







    @staticmethod
    def filter_low_pwr_lags(raw_data,range_list):
        """Removes low power lags from fitting

        Prunes off low power lags determined by cutoff criteria. Once a
        cutoff lag is determined, all subsequent lags in the list are 
        removed

        :param raw_data: a dictionary of raw data parameters
        :param range_list: A list of Range objects with data points

        """
        pwr0 = raw_data['pwr0']
        mplgs = raw_data['mplgs']
        nave = raw_data['nave']

        for range_obj in range_list:
            range_number = range_obj.range_number
            log_sigma_fluc = math.log(FLUCTUATION_CUTOFF_COEFF * pwr0[range_number]/math.sqrt(2 * nave))

            bad_indices = []
            cutoff_lag = mplgs + 1

            # iterator = np.nditer([range_obj.pwrs.log_pwrs,range_obj.alpha_2],flags=['f_index'])
            # while not iterator.finished:
            #   log_pwr,alpha_2 = iterator[0],iterator[1]
            #   idx = iterator.index
            #   if idx > cutoff_lag:
            #       bad_indices.append(idx)
            #   else:
            #       if((1/math.sqrt(alpha_2) <= ALPHA_CUTOFF) and (log_pwr <= log_sigma_fluc)):
            #           cutoff_lag = idx
            #           bad_indices.append(idx)

            #   iterator.iternext()

            for idx,(log_pwr,alpha_2) in enumerate(zip(range_obj.pwrs.log_pwrs,range_obj.alpha_2)):
                if idx > cutoff_lag:
                    bad_indices.append(idx)
                else:
                    if((1/math.sqrt(alpha_2) <= ALPHA_CUTOFF) and (log_pwr <= log_sigma_fluc)):
                        cutoff_lag = idx
                        bad_indices.append(idx)

            range_obj.pwrs.remove_bad_points(bad_indices)


    @staticmethod
    def filter_bad_acfs(raw_data,ranges,noise_pwr):
        """Removes bad ACFs entirely from analysis

        Removes ACFs which are deemed to be pure noise, or ACFs with too
        few power lags left

        :param raw_data: a dictionary of raw data parameters
        :param range_list: A list of Range objects with data points
        :param noise_pwr: minimum power for which an ACF is pure noise
        
        """
        nave = raw_data['nave']
        pwr0 = raw_data['pwr0']

        cutoff_pwr = ACF_SNR_CUTOFF * noise_pwr * (1 + 1/math.sqrt(nave))

        bad_indices = []
        for idx,rang in enumerate(ranges):
            range_number = rang.range_number
            if (pwr0[range_number] <= cutoff_pwr) or (len(rang.pwrs.log_pwrs) < MIN_LAGS):
                bad_indices.append(idx)

        for i in sorted(bad_indices,reverse=True):
            del ranges[i]

    

class Determinations(object):
    """This is a class to construct a new dictionary of final deteminations

    This class holds the methods to convert fitted data into the final measurements
    of the plasma

    """

    def __init__(self,raw_data,range_list,noise_pwr):
        hdw_info = hdw[raw_data['stid']]
        self.paramater_dict = self.new_parameter_dictionary(hdw_info,raw_data,range_list,noise_pwr)
    
    def new_parameter_dictionary(self,hdw_info,raw_data,range_list,noise_pwr):
        """Creates the new dictionary of parameters from fitted data

        :param hdw_info: a dictionary of the radar hardware information
        :param raw_data: a dictionary of raw data parameters
        :param range_list: a list of Range objects with fitted data
        :param noise_pwr: minimum power for which an ACF is pure noise
        :returns: a new dictionary of the fitacf parameters

        """

        new_parameter_dict = {}

        number_of_ranges = len(range_list)

        new_parameter_dict['radar.revision.major'] = raw_data['radar.revision.major']
        new_parameter_dict['radar.revision.minor'] = raw_data['radar.revision.minor']
        new_parameter_dict['origin.code'] = raw_data['origin.code']
        new_parameter_dict['origin.time'] = datetime.utcnow().strftime('%c')
        new_parameter_dict['origin.command'] = " ".join(sys.argv)
        new_parameter_dict['cp'] = raw_data['cp']
        new_parameter_dict['stid'] = raw_data['stid']
        new_parameter_dict['time.yr'] = raw_data['time.yr']
        new_parameter_dict['time.mo'] = raw_data['time.mo']
        new_parameter_dict['time.dy'] = raw_data['time.dy']
        new_parameter_dict['time.hr'] = raw_data['time.hr']
        new_parameter_dict['time.mt'] = raw_data['time.mt']
        new_parameter_dict['time.sc'] = raw_data['time.sc']
        new_parameter_dict['time.us'] = raw_data['time.us']
        new_parameter_dict['txpow'] = raw_data['txpow']
        new_parameter_dict['nave'] = raw_data['nave']
        new_parameter_dict['atten'] = raw_data['atten']
        new_parameter_dict['lagfr'] = raw_data['lagfr']
        new_parameter_dict['smsep'] = raw_data['smsep']
        new_parameter_dict['ercod'] = raw_data['ercod']
        new_parameter_dict['stat.agc'] = raw_data['stat.agc']
        new_parameter_dict['stat.lopwr'] = raw_data['stat.lopwr']
        new_parameter_dict['noise.search'] = raw_data['noise.search']
        new_parameter_dict['noise.mean'] = raw_data['noise.mean']
        new_parameter_dict['channel'] = raw_data['channel']
        new_parameter_dict['bmnum'] = raw_data['bmnum']
        new_parameter_dict['bmazm'] = raw_data['bmazm']
        new_parameter_dict['scan'] = raw_data['scan']
        new_parameter_dict['offset'] = raw_data['offset']
        new_parameter_dict['rxrise'] = raw_data['rxrise']
        new_parameter_dict['intt.sc'] = raw_data['intt.sc']
        new_parameter_dict['intt.us'] = raw_data['intt.us']
        new_parameter_dict['txpl'] = raw_data['txpl']
        new_parameter_dict['mpinc'] = raw_data['mpinc']
        new_parameter_dict['mppul'] = raw_data['mppul']
        new_parameter_dict['mplgs'] = raw_data['mplgs']
        new_parameter_dict['nrang'] = raw_data['nrang']
        new_parameter_dict['frang'] = raw_data['frang']
        new_parameter_dict['rsep'] = raw_data['rsep']
        new_parameter_dict['xcf'] = raw_data['xcf']
        new_parameter_dict['tfreq'] = raw_data['tfreq']
        new_parameter_dict['mxpwr'] = raw_data['mxpwr']
        new_parameter_dict['lvmax'] = raw_data['lvmax']
        new_parameter_dict['fitacf.revision.major'] = FITACF_REVISION_MAJOR
        new_parameter_dict['fitacf.revision.minor'] = FITACF_REVISON_MINOR
        new_parameter_dict['combf'] = raw_data['combf']
        new_parameter_dict['noise.sky'] = noise_pwr
        new_parameter_dict['noise.lag0'] = 0.0
        new_parameter_dict['noise.vel'] = 0.0
        new_parameter_dict['ptab'] = raw_data['ptab']
        new_parameter_dict['ltab'] = raw_data['ltab']
        new_parameter_dict['pwr0'] = self.lag_0_pwr_in_dB(raw_data,noise_pwr)
        new_parameter_dict['slist'] = self.set_slist(range_list)
        new_parameter_dict['nlag'] = self.set_nlag(range_list)
        new_parameter_dict['qflg'] = self.set_qflg(range_list)
        new_parameter_dict['p_l'] = self.set_p_l(range_list,noise_pwr)
        new_parameter_dict['p_l_e'] = self.set_p_l_err(range_list)
        new_parameter_dict['p_s'] = self.set_p_s(range_list,noise_pwr)
        new_parameter_dict['p_s_e'] = self.set_p_s_err(range_list)
        new_parameter_dict['v'] = self.set_v(range_list,raw_data,hdw_info)
        new_parameter_dict['v_e'] = self.set_v_err(range_list,raw_data,hdw_info)
        new_parameter_dict['w_l'] = self.set_w_l(range_list,raw_data)
        new_parameter_dict['w_l_e'] = self.set_w_l_err(range_list,raw_data)
        new_parameter_dict['w_s'] = self.set_w_s(range_list,raw_data)
        new_parameter_dict['w_s_e'] = self.set_w_s_err(range_list,raw_data)
        new_parameter_dict['sd_l'] = self.set_sdev_l(range_list)
        new_parameter_dict['sd_s'] = self.set_sdev_s(range_list)
        new_parameter_dict['sd_phi'] = self.set_sdev_phi(range_list)
        new_parameter_dict['gflg'] = self.set_gsct(new_parameter_dict['v'],new_parameter_dict['w_l'])
        new_parameter_dict['x_qflg'] = [0] * number_of_ranges
        new_parameter_dict['x_gflg'] = [0] * number_of_ranges
        new_parameter_dict['x_p_l'] = [0] * number_of_ranges
        new_parameter_dict['x_p_l_e'] = [0] * number_of_ranges
        new_parameter_dict['x_p_s'] = [0] * number_of_ranges
        new_parameter_dict['x_p_s_e'] = [0] * number_of_ranges
        new_parameter_dict['x_v'] = [0] * number_of_ranges
        new_parameter_dict['x_v_e'] = [0] * number_of_ranges
        new_parameter_dict['x_w_l'] = [0] * number_of_ranges
        new_parameter_dict['x_w_l_e'] = [0] * number_of_ranges
        new_parameter_dict['x_w_s'] = [0] * number_of_ranges
        new_parameter_dict['x_w_s_e'] = [0] * number_of_ranges
        new_parameter_dict['phi0'] = self.set_xcf_phi0(range_list)
        new_parameter_dict['phi0_e'] = self.set_xcf_phi0_err(range_list)

        elv = self.find_elevation(range_list,raw_data,hdw_info)
        new_parameter_dict['elv_low'] = elv[0]
        new_parameter_dict['elv'] = elv[1]
        new_parameter_dict['elv_high'] = elv[2]
        new_parameter_dict['x_sd_l'] = [0] * number_of_ranges
        new_parameter_dict['x_sd_s'] = [0] * number_of_ranges
        new_parameter_dict['x_sd_phi'] = self.set_xcf_sdev_phi(range_list)

        return new_parameter_dict


    def lag_0_pwr_in_dB(self,raw_data,noise_pwr):
        """Converts lag 0 powers to dB

        :param raw_data: a dictionary of raw data parameters
        :param noise_pwr: minimum power for which an ACF is pure noise
        :returns: an array of lag 0 powers in dB

        """

        pwr0 = raw_data['pwr0']

        pwr_conversion = lambda x: 10 * np.log10((x - noise_pwr)/noise_pwr)
        pwr_dB = [pwr_conversion(pwr) if (pwr - noise_pwr > 0.0) else -50.0 for pwr in pwr0]

        return np.array(pwr_dB)


    def set_slist(self,range_list):
        """Creates the array of good ranges

        :param range_list: a list of Range objects left after filtering
        :returns: an array of range numbers left after filtering

        """

        slist = [range_obj.range_number for range_obj in range_list]

        return np.array(slist)


    def set_nlag(self,range_list):
        """Sets the number of points used for power fitting

        :param range_list: a list of Range objects with data points
        :returns: an array of the number of points left for fitting at each range

        """

        nlag = [len(range_obj.pwrs.log_pwrs) for range_obj in range_list]

        return np.array(nlag)


    def set_qflg(self,range_list):
        """Creates a qflg array

        All data is valid at this point so this
        just makes an array of ones the length of the range list.

        :param range_list: a list of Range objects after filtering
        :returns: an array of ones
        """

        return np.ones(len(range_list),dtype=np.int64)


    def set_p_l(self,range_list,noise_pwr):
        """Computes the power in dB of the linear power fit at each range

        :param range_list: a list of Range objects after fitting
        :param noise_pwr: minimum power for which an ACF is pure noise
        :returns: an array of fitted lambda powers in dB

        """

        noise_dB = 10 * math.log10(noise_pwr);

        p_l_conversion = lambda x: 10 * x/math.log(10) - noise_dB

        p_l = [p_l_conversion(range_obj.linear_pwr_fit.a) for range_obj in range_list]

        return np.array(p_l)


    def set_p_l_err(self,range_list):
        """Computes the power in dB of the linear power fit error at each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted lambda power errors in dB

        """

        p_l_err_conversion = lambda x: 10 * math.sqrt(x) / math.log(10)

        p_l_err = [p_l_err_conversion(range_obj.linear_pwr_fit_err.sigma_2_a) for range_obj in range_list]

        return np.array(p_l_err)


    def set_p_s(self,range_list,noise_pwr):
        """Computes the power in dB of the quadratic power fit at each range

        :param range_list: a list of Range objects after fitting
        :param noise_pwr: minimum power for which an ACF is pure noise
        :returns: an array of fitted sigma powers in dB

        """

        noise_dB = 10 * math.log10(noise_pwr)

        p_s_conversion = lambda x: 10 * x/math.log(10) - noise_dB

        p_s = [p_s_conversion(range_obj.quadratic_pwr_fit.a) for range_obj in range_list]

        return np.array(p_s)


    def set_p_s_err(self,range_list):
        """Computes the power in dB of the quadratic power fit error at each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted sigma power errors in dB

        """

        p_s_err_conversion = lambda x: 10 * math.sqrt(x) / math.log(10)

        p_s_err = [p_s_err_conversion(range_obj.quadratic_pwr_fit_err.sigma_2_a) for range_obj in range_list]

        return np.array(p_s_err)


    def set_v(self,range_list,raw_data,hdw_info):
        """Computes the fitted velocity at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :param hdw_info: a dictionary of the radar hardware information
        :returns: an array of computed velocites

        """

        vel_conversion = C/((4*math.pi)*(raw_data['tfreq'] * 1000.0)) * hdw_info['velsign']

        vel_calculation = lambda x: x * vel_conversion
        vel = [vel_calculation(range_obj.phase_fit.b) for range_obj in range_list]

        return np.array(vel)


    def set_v_err(self,range_list,raw_data,hdw_info):
        """Computes the fitted velocity error at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :param hdw_info: a dictionary of the radar hardware information
        :returns: an array of computed velocites

        """        
        vel_conversion = C/((4*math.pi)*(raw_data['tfreq'] * 1000.0)) * hdw_info['velsign']

        vel_err = [math.sqrt(range_obj.phase_fit.sigma_2_b) * vel_conversion for range_obj in range_list]

        return np.array(vel_err)


    def set_w_l(self,range_list,raw_data):
        """Computes the spectral width from the linear power fit at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral widths

        """

        width_conversion = C/((4*math.pi)*(raw_data['tfreq'] * 1000.0))*2.0

        w_l_calculation = lambda x: math.fabs(x) * width_conversion
        w_l = [w_l_calculation(range_obj.linear_pwr_fit.b) for range_obj in range_list]

        return np.array(w_l)
        

    def set_w_l_err(self,range_list,raw_data):
        """Computes the spectral width error from the linear power fit errors
        at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral width errors

        """

        width_conversion = C/((4*math.pi)*(raw_data['tfreq'] * 1000.0))*2.0

        w_l_err_calculation = lambda x: math.sqrt(x) * width_conversion
        w_l_err = [w_l_err_calculation(range_obj.linear_pwr_fit_err.sigma_2_b) for range_obj in range_list]

        return np.array(w_l_err)

    def set_w_s(self,range_list,raw_data):
        """Computes the spectral width from the quadratic power fit at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral widths

        """        

        w_s_conversion = C/(4*math.pi)/(raw_data['tfreq'] * 1000.0) *4.* math.sqrt(math.log(2))

        w_s_calculation = lambda x: math.sqrt(math.fabs(x)) * w_s_conversion
        w_s = [w_s_calculation(range_obj.quadratic_pwr_fit.b) for range_obj in range_list]

        return np.array(w_s)

    def set_w_s_err(self,range_list,raw_data):
        """Computes the spectral width error from the quadratic power fit errors
        at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral width errors

        """

        w_s_conversion = C/(4*math.pi)/(raw_data['tfreq'] * 1000.0) *4.* math.sqrt(math.log(2))

        w_s_calculation = lambda x,y: math.sqrt(math.fabs(x))/2./math.sqrt(math.fabs(y)) * w_s_conversion
        w_s_err = [w_s_calculation(range_obj.quadratic_pwr_fit_err.sigma_2_b,range_obj.quadratic_pwr_fit.b) for range_obj in range_list]        

        return np.array(w_s_err)


    def set_sdev_l(self,range_list):
        """Sets the chi_2 value of the linear power fit for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of chi_2 values

        """

        s_sdev_l = [range_obj.linear_pwr_fit.chi_2 for range_obj in range_list]

        return np.array(s_sdev_l)


    def set_sdev_s(self,range_list):
        """Sets the chi_2 value of the quadratic power fit for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of chi_2 values

        """

        s_sdev_s = [range_obj.quadratic_pwr_fit.chi_2 for range_obj in range_list]

        return np.array(s_sdev_s)


    def set_sdev_phi(self,range_list):
        """Sets the chi_2 value of the ACF phase fit for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of chi_2 values

        """
        s_sdev_phi = [range_obj.phase_fit.chi_2 for range_obj in range_list]

        return np.array(s_sdev_phi)

    def set_gsct(self,vel_arr,w_l_arr):
        """Computes whether scatter comes from ground or ionsphere for each
        range

        :param vel_arr: an array of determined velocities
        :param w_l_arr: an array of determined spectral widths(lambda)
        :returns: an array of ground scatter flags

        """
        v_abs = np.fabs(vel_arr)
        gsct_criteria = v_abs - (V_MAX - w_l_arr * (V_MAX/W_MAX))

        gsct = [1 if g < 0.0 else 0 for g in gsct_criteria]

        return np.array(gsct)

    def find_elevation(self,range_list,raw_data,hdw_info):
        """Computes elevation angle for each range

        Computes fitted elevation angle, unfitted elevation angle, and
        error in elevation all in one method

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :param hdw_info: a dictionary of the radar hardware information
        :returns: a tuple of arrays with errors,fitted elevation, and unfitted elevation

        """
        bmnum = raw_data['bmnum']
        tfreq = raw_data['tfreq']

        x = hdw_info['interoffx']
        y = hdw_info['interoffy']
        z = hdw_info['interoffz']

        antenna_sep = math.sqrt(x*x + y*y + z*z)

        elev_corr = hdw_info['phasesign'] * math.asin(z/antenna_sep)

        if y > 0.0:
            phi_sign = 1
        else:
            phi_sign = -1
            elev_corr = -elev_corr

        azi_offset = hdw_info['maxbeams']/2 - 0.5
        phi_0 = hdw_info['beamsep'] * ( bmnum- azi_offset) * math.pi/180
        c_phi_0 = math.cos(phi_0)

        wave_num = 2 * math.pi * tfreq * 1000/C;

        cable_offset = -2 * math.pi * tfreq * 1000 * hdw_info['tdiff'] * 1.0e-6

        phase_diff_max = phi_sign * wave_num * antenna_sep * c_phi_0 + cable_offset

        psi_calculation = lambda x: x + 2 * math.pi * math.floor((phase_diff_max-x)/(2*math.pi))
        psi_uncorrected = [psi_calculation(range_obj.elev_fit.a) for range_obj in range_list]

        if(phi_sign < 0):
            psi_uncorrected = [psi_u + 2 * math.pi for psi_u in psi_uncorrected]
            #psi_uncorrected = psi_uncorrected + 2 * math.pi
        #psi_uncorrected = np.array(psi_uncorrected)

        psi = [psi_u - cable_offset for psi_u in psi_uncorrected]

        psi_kd = [p/(wave_num * antenna_sep) for p in psi]
        theta = [c_phi_0**2 - pkd**2 for pkd in psi_kd] 

        # if( (theta < 0.0) or (math.fabs(theta) > 1.0)):
        #   elevation = -elev_corr
        # else:
        #   elevation = math.asin(math.sqrt(theta))

        elev_calculation = lambda x: -elev_corr if (t < 0.0 or math.fabs(x) > 1.0) else math.asin(math.sqrt(x)) 
        elevation = [elev_calculation(t) for t in theta]
        
        elevation_normal = [180/math.pi * (elev + elev_corr) for elev in elevation]

        #Elevation errors
        psi_k2d2 = [p/(wave_num**2 * antenna_sep**2) for p in psi]
        df_by_dy = [pkd/math.sqrt(t * (1 - t)) for pkd,t in zip(psi_k2d2,theta)]
        
        elev_low_calculation = lambda x,y: 180/math.pi * math.sqrt(x) * math.fabs(y)
        elevation_low = [elev_low_calculation(range_obj.elev_fit.sigma_2_a,dfdy) for range_obj,dfdy in zip(range_list,df_by_dy)]


        #Experiment to compare fitted and measured elevation
        xcfd = raw_data['xcfd']
        real = [xcfd[range_obj.range_number][0][0] for range_obj in range_list]
        imag = [xcfd[range_obj.range_number][0][1] for range_obj in range_list]
        xcf0_p = [math.atan2(i,r) for i,r in zip(imag,real)]

        psi_uu_calculation = lambda x: x + 2 * math.pi * math.floor((phase_diff_max-x)/(2*math.pi))
        psi_uncorrected_unfitted = [psi_uu_calculation(x) for x in xcf0_p]

        # if phi_sign < 0:
        #   psi_uncorrected_unfitted = psi_uncorrected_unfitted + (2 * math.pi)

        psi_uu_calculation = lambda x: x + (2 * math.pi) if phi_sign < 0 else x
        psi_uncorrected_unfitted = [psi_uu_calculation(p_uu) for p_uu in psi_uncorrected_unfitted]
        
        psi = [p_uu - cable_offset for p_uu in psi_uncorrected_unfitted]

        psi_kd = [p/(wave_num * antenna_sep) for p in psi]
        theta = [c_phi_0**2 - pkd**2 for pkd in psi_kd]

        # if( (theta < 0.0) or (np.fabs(theta) > 1.0) ){
        #   elevation = -elev_corr
        # }
        # else{
        #   elevation = np.asin(np.sqrt(theta))
        # }

        elev_calculation = lambda x: -elev_corr if (t < 0.0 or math.fabs(x) > 1.0) else math.asin(math.sqrt(x)) 
        elevation = [elev_calculation(t) for t in theta]

        elevation_high = [180/math.pi * (elev + elev_corr) for elev in elevation]

        return (elevation_low,elevation_normal,elevation_high)


    def set_xcf_phi0(self,range_list):
        """Sets the fitted offset of the XCF phase for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted XCF phases offsets

        """
        phi0 = [range_obj.elev_fit.a for range_obj in range_list]

        return np.array(phi0)


    def set_xcf_phi0_err(self,range_list):
        """Sets the fitted offset error of the XCF phase for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted XCF phases offset errors

        """
        phi0_err = [math.sqrt(range_obj.elev_fit.sigma_2_a) for range_obj in range_list]

        return np.array(phi0_err)

    def set_xcf_sdev_phi(self,range_list):
        """Sets the chi_2 value for the XCF phase fit for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of chi_2 values

        """
        sdev_phi = [range_obj.elev_fit.chi_2 for range_obj in range_list]

        return np.array(sdev_phi) 



def _fit(raw_data):
    """Fits a single dictionary of raw data.

    This function is meant to be passed as an argument
    to multiprocessing map.

    :param raw_data: a dictionary of raw data parameters
    :returns: a dictionary of fitted parameters

    """

    lags = create_lag_list(raw_data)
    noise_pwr = NoisePower.acf_cutoff_pwr(raw_data)

    range_list = []
    for i in range(0,raw_data['nrang']):
        if raw_data['pwr0'][i] != 0:
            new_range = Range(i,raw_data,lags)
            range_list.append(new_range)


    Filtering.filter_tx_overlapped_lags(raw_data,lags,range_list)
    Filtering.filter_inf_lags(range_list)
    Filtering.filter_low_pwr_lags(raw_data,range_list)
    Filtering.filter_bad_acfs(raw_data,range_list,noise_pwr)

    ACFFitting.acf_pwr_fitting(range_list)

    ACFFitting.calculate_phase_and_elev_sigmas(range_list,raw_data)

    ACFFitting.acf_phase_unwrap(range_list)
    ACFFitting.acf_phase_fitting(range_list)

    ACFFitting.xcf_phase_unwrap(range_list)
    ACFFitting.xcf_phase_fitting(range_list)

    determined_parameters = Determinations(raw_data,range_list,noise_pwr)
    return determined_parameters.paramater_dict


def fit(raw_records):
    """Performs the whole fitting procedure for rawacf data

    Calls the _fit procedure in a parallelized multiprocessing
    environment to speed up the procedure. The speed of this
    routine scales with number of cores.

    :param raw_records: a list of raw data dictionaries
    :returns: a list of dictionaries with fitted data

    """

    pool = mp.Pool()
    fitted_records = pool.map(_fit,raw_records)

    return fitted_records




if __name__ == "__main__":
    in_filename = sys.argv[1]
    raw_records = dm.parse_dmap_format_from_file(in_filename)
    
    fitted_records = fit(raw_records)

    # for rr in raw_records:
    #     _fit(rr)

    out_filename = sys.argv[2]
    dm.dicts_to_file(fitted_records,out_filename,file_type='fitacf')