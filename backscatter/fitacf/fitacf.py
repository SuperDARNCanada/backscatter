from __future__ import print_function
import numpy as np
import sys
import math
#import mpmath as mpm
from datetime import datetime
import multiprocessing as mp
import warnings
import argparse
from backscatter import dmap as dm

from backscatter.fitacf.noisepower import NoisePower
from backscatter.fitacf.filtering import Filtering
from backscatter.fitacf.fitting import ACFFitting
from backscatter.fitacf.determinations import Determinations


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
        self.alpha_2 = range_obj.alpha_2
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
        acfd = raw_data['acfd'][range_obj.range_idx]
        mplgs = raw_data['mplgs']
        pwr0 = raw_data['pwr0'][range_obj.range_number]
        nave = raw_data['nave']
        mpinc = raw_data['mpinc']

        real = acfd[:,0]
        imag = acfd[:,1]
        real_2 = real**2
        imag_2 = imag**2

        pwrs = np.sqrt(np.add(real_2,imag_2))

        pwr_normalized = pwrs/pwr0
        pwr_normalized_2 = pwr_normalized**2
        inverse_alpha_2 = np.reciprocal(range_obj.alpha_2)

        sigmas = [pwr0 * np.sqrt((pn_2 + ia_2)/(2 * nave)) for (pwr, pn_2, ia_2) in zip(pwrs, pwr_normalized_2,inverse_alpha_2)]

        self.sigmas = np.array(sigmas)
        self.t = np.round(np.array([lag['number'] * mpinc * 1.0e-6 for lag in lags]),decimals=6)

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
        self.alpha_2 = self.alpha_2[mask]

    def remove_inf_points(self,non_inf_indices):

        self.log_pwrs = self.log_pwrs[non_inf_indices]
        self.sigmas = self.sigmas[non_inf_indices]
        self.t = self.t[non_inf_indices]
        self.alpha_2 = self.alpha_2[mask]


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
        self.alpha_2 = range_obj.alpha_2

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

        mplgs = raw_data['mplgs']
        mpinc = raw_data['mpinc']

        if phase_type not in raw_data:
            acfd = np.zeros((mplgs,2))
        else:
            acfd = raw_data[phase_type][range_obj.range_idx]

        real = acfd[:,0]
        imag = acfd[:,1]

        self.phases = np.arctan2(imag,real)

        self.sigmas = np.zeros(mplgs)

        self.t = np.round(np.array([lag['number'] * mpinc * 1.0e-6 for lag in lags]),decimals=6)

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
        self.alpha_2 = self.alpha_2[mask]

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

    def __init__(self,idx,range_number,raw_data,lags):
        self.range_number = range_number
        self.range_idx = idx
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
                range_to_check = int(pulse_diff * tau + self.range_number)

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

def debug_output(range_list, raw_data):
    for range_obj in range_list:
        time = "TIME %d-%02d-%02dT%02d:%02d:%f"%(raw_data['time.yr'], raw_data['time.mo'],
                          raw_data['time.dy'], raw_data['time.hr'],
                          raw_data['time.mt'], raw_data['time.sc'] +
                          raw_data['time.us']/1.0e6)
        print(time)
        print("RANGE NUM", range_obj.range_number)
        print("linear_pwr_fit.a", range_obj.linear_pwr_fit.a)
        print("linear_pwr_fit.b", range_obj.linear_pwr_fit.b)
        print("linear_pwr_fit.sigma_2_a",range_obj.linear_pwr_fit.sigma_2_a)
        print("linear_pwr_fit.sigma_2_b",range_obj.linear_pwr_fit.sigma_2_b)
        print("quadratic_pwr_fit.a",range_obj.quadratic_pwr_fit.a)
        print("quadratic_pwr_fit.b",range_obj.quadratic_pwr_fit.b)
        print("quadratic_pwr_fit.sigma_2_a",range_obj.quadratic_pwr_fit.sigma_2_a)
        print("quadratic_pwr_fit.sigma_2_b",range_obj.quadratic_pwr_fit.sigma_2_b)
        print("linear_pwr_fit_err.a",range_obj.linear_pwr_fit_err.a)
        print("linear_pwr_fit_err.b",range_obj.linear_pwr_fit_err.b)
        print("linear_pwr_fit_err.sigma_2_a",range_obj.linear_pwr_fit_err.sigma_2_a)
        print("linear_pwr_fit_err.sigma_2_b",range_obj.linear_pwr_fit_err.sigma_2_b)
        print("quadratic_pwr_fit_err.a",range_obj.quadratic_pwr_fit_err.a)
        print("quadratic_pwr_fit_err.b",range_obj.quadratic_pwr_fit_err.b)
        print("quadratic_pwr_fit_err.sigma_2_a",range_obj.quadratic_pwr_fit_err.sigma_2_a)
        print("quadratic_pwr_fit_err.sigma_2_b",range_obj.quadratic_pwr_fit_err.sigma_2_b)
        print("alphas", range_obj.alpha_2)
        print("ACF LOG POWERS")
        for p,s,t in zip(range_obj.pwrs.log_pwrs, range_obj.pwrs.sigmas,range_obj.pwrs.t):
            print("power",p,"sigma",s,"t",t)
        print("ACF PHASES")
        for p,s,t in zip(range_obj.phases.phases,range_obj.phases.sigmas,range_obj.phases.t):
            print("phi",p,"sigma",s,"t",t)
        print("phase_fit.a",range_obj.phase_fit.a)
        print("phase_fit.b",range_obj.phase_fit.b)
        print("phase_fit.sigma_2_a",range_obj.phase_fit.sigma_2_a)
        print("phase_fit.sigma_2_b",range_obj.phase_fit.sigma_2_b)
        print("XCF PHASES")
        for p,s,t in zip(range_obj.elevs.phases,range_obj.elevs.sigmas,range_obj.elevs.t):
            print("phi",p,"sigma",s,"t",t)
        print("elev_fit.S", range_obj.elev_fit.S)
        print("elev_fit.S_x", range_obj.elev_fit.S_x)
        print("elev_fit.S_y", range_obj.elev_fit.S_y)
        print("elev_fit.S_xx", range_obj.elev_fit.S_xx)
        print("elev_fit.S_xy", range_obj.elev_fit.S_xy)
        print("elev_fit.delta", range_obj.elev_fit.delta)
        print("elev_fit.delta_a", range_obj.elev_fit.delta_a)
        print("elev_fit.delta_b", range_obj.elev_fit.delta_b)
        print("elev_fit.a",range_obj.elev_fit.a)
        print("elev_fit.b",range_obj.elev_fit.b)
        print("elev_fit.sigma_2_a",range_obj.elev_fit.sigma_2_a)
        print("elev_fit.sigma_2_b",range_obj.elev_fit.sigma_2_b)
        print("elev_fit.cov_ab", range_obj.elev_fit.cov_ab)
        print("elev_fit.r_ab", range_obj.elev_fit.r_ab)
        print("elev_fit.Q", range_obj.elev_fit.Q)
        print("elev_fit.chi_2", range_obj.elev_fit.chi_2)

def _fit(raw_data, debug_mode=False):
    """Fits a single dictionary of raw data.

    This function is meant to be passed as an argument
    to multiprocessing map.

    :param raw_data: a dictionary of raw data parameters
    :returns: a dictionary of fitted parameters

    """

    lags = create_lag_list(raw_data)

    # We check number of averages < 0 since this will cause invalid
    # division in the noise calculation
    if raw_data['nave'] <= 0:
        noise_pwr = 1.0
    else:
        raw_data['nave']
        noise_pwr = NoisePower.acf_cutoff_pwr(raw_data)

    range_list = []
    #num_ranges_with_data = len(raw_data['slist'])
    for idx,range_number in enumerate(raw_data['slist']):
        if raw_data['pwr0'][range_number] != 0:
            new_range = Range(idx,range_number,raw_data,lags)
            range_list.append(new_range)


    Filtering.filter_tx_overlapped_lags(raw_data,lags,range_list)
    Filtering.filter_inf_lags(range_list)
    Filtering.filter_low_pwr_lags(raw_data,range_list)
    Filtering.filter_bad_acfs(raw_data,range_list,noise_pwr)

    ACFFitting.acf_pwr_fitting(range_list)

    ACFFitting.calculate_phase_and_elev_sigmas(range_list,raw_data)

    ACFFitting.acf_phase_unwrap(range_list,raw_data)
    ACFFitting.acf_phase_fitting(range_list)

    Filtering.filter_bad_fits(range_list)

    ACFFitting.xcf_phase_unwrap(range_list)
    ACFFitting.xcf_phase_fitting(range_list)

    if debug_mode:
        debug_output(range_list, raw_data)

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
    parser = argparse.ArgumentParser()
    parser.add_argument("infile", help="Input rawacf file path", type=str)
    parser.add_argument("outfile", help="Output file path", type=str)
    parser.add_argument("--debug-mode", help="""Disables multiprocessing usage to produce more
        meaningful exception at the cost of performance.""", action="store_true")
    args = parser.parse_args()

    raw_records = dm.parse_dmap_format_from_file(args.infile)

    if args.debug_mode:
        fitted_records = [_fit(rr,True) for rr in raw_records]
    else:
        fitted_records = fit(raw_records)

    dm.dicts_to_file(fitted_records,args.outfile,file_type='fitacf')
