from backscatter import config
import numpy as np
import math

MIN_LAGS = int(config.get("fitacf","minimum_lags"))
ACF_SNR_CUTOFF = float(config.get("fitacf","acf_snr_cutoff"))
ALPHA_CUTOFF = float(config.get("fitacf","alpha_cutoff"))
FLUCTUATION_CUTOFF_COEFF = int(config.get("fitacf","fluctuation_cutoff_coeff"))

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
        offset = raw_data['offset']
        channel = raw_data['channel']

        pulses_in_us = [mpinc * pulse for pulse in pulses]

        pulses_stereo = []
        if offset != 0 and (channel == 1 or channel == 2):
            for pulse in pulses:
                if channel == 1:
                    pulse_us = pulse * mpinc - offset
                else:
                    pulse_us = pulse * mpinc + offset

                pulses_stereo.append(pulse_us)

        pulses_in_us = pulses_in_us + pulses_stereo
        i = -1
        ts, t1, t2, sample = lagfr, 0, 0, 0
        bad_samples = []
        for pulse_us in pulses_in_us:

            t1 = pulse_us - txpl/2
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

        #Division by zero error
        if nave <= 0:
            return

        for range_obj in range_list:
            range_number = range_obj.range_number
            if len(range_obj.pwrs.log_pwrs) == 0:
                continue

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
            noise_mean = raw_data['noise.mean']
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

        #Division by zero error
        if nave <= 0:
            return

        cutoff_pwr = noise_pwr * 2

        bad_indices = []
        for idx,rang in enumerate(ranges):
            range_number = rang.range_number
            pw = float(pwr0[range_number])
            ln = len(rang.pwrs.log_pwrs)
            if (pw <= cutoff_pwr) or (ln < MIN_LAGS):
                bad_indices.append(idx)


        for i in sorted(bad_indices,reverse=True):
            del ranges[i]
