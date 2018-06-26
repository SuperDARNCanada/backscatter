from backscatter import config
import numpy as np
import math

MIN_LAGS = int(config.get("fitacf", "minimum_lags"))  # REVIEW #7 Documentation for these globals?
ACF_SNR_CUTOFF = float(config.get("fitacf", "acf_snr_cutoff"))   # REVIEW #7 What is the ACR SNR Cutoff?
ALPHA_CUTOFF = float(config.get("fitacf", "alpha_cutoff"))  # REVIEW #7 What is the alpha cutoff?
FLUCTUATION_CUTOFF_COEFF = int(config.get("fitacf", "fluctuation_cutoff_coeff"))  # REVIEW #7 What is the fluctuation cutoff coefficient?


class Filtering(object):
    """This class contains static methods relating to the filtering
    of bad data from the consideration of fitting.

    """

    @staticmethod
    def mark_bad_samples(raw_data):
        """Mark the samples that are blacked out by TX overlap

        :param raw_data: a dictionary of raw data parameters

        """
        mppul = raw_data['mppul']  # REVIEW #42 - mppul is not used
        pulses = raw_data['ptab']
        mpinc = raw_data['mpinc']
        txpl = raw_data['txpl']
        smsep = raw_data['smsep']
        lagfr = raw_data['lagfr']
        offset = raw_data['offset']
        channel = raw_data['channel']

        pulses_in_us = [mpinc * pulse for pulse in pulses]

        pulses_stereo = []  # REVIEW #2 Maybe a comment about twofsound and how it deals with channels... this seems like it will handle the twofsound files correctly since the if statment won't be entered.
        if offset != 0 and (channel == 1 or channel == 2):  # REVIEW #40 - consistency with parentheses
            for pulse in pulses:
                if channel == 1:
                    pulse_us = pulse * mpinc - offset
                else:
                    pulse_us = pulse * mpinc + offset

                pulses_stereo.append(pulse_us)

        pulses_in_us = pulses_in_us + pulses_stereo
        i = -1  # REVIEW #42 - 'i' is not used
        ts, t1, t2, sample = lagfr, 0, 0, 0
        bad_samples = []
        for pulse_us in pulses_in_us:
            # REVIEW #29 - what are 2, 3 and 100 below?
            t1 = pulse_us - txpl / 2 # REVIEW #3- This looks like it assumes that the samples are occurring right in the middle of the pulses, is that the case?
            t2 = t1 + 3 * txpl / 2 + 100

            # we now have a pulse that occurs after the current sample.  Start
            # incrementing the sample number until we find a sample that lies
            # within the pulse
            while (ts < t1):  # REVIEW # - consistency with parentheses
                sample = sample + 1
                ts = ts + smsep

            # ok, we now have a sample which occurs after the pulse starts.
            # check to see if it occurs before the pulse ends, and if so, mark
            # it as a bad sample
            while ((ts >= t1) and (ts <= t2)):  # REVIEW #40 - consistency with parentheses
                bad_samples.append(sample)
                sample = sample + 1
                ts = ts + smsep

        return bad_samples

    @staticmethod
    def filter_tx_overlapped_lags(raw_data, lags, range_list):
        """Remove data points affected by TX overlapped lags

        :param raw_data: a dictionary of raw data parameters
        :param lags: list of lag dictionaries
        :param range_list: A list of Range objects with data points

        """

        bad_samples = Filtering.mark_bad_samples(raw_data)

        for range_obj in range_list:

            bad_indices = []
            for idx, lag in enumerate(lags):
                sample1 = lag['sample_base1'] + range_obj.range_number # REVIEW #26 Are 'sample_base1' and 'sample_base2' good names? We're not sure, starting with 0 would be better as well
                sample2 = lag['sample_base2'] + range_obj.range_number

                if (sample1 in bad_samples) or (sample2 in bad_samples):
                    bad_indices.append(idx)

            range_obj.pwrs.remove_bad_points(bad_indices)
            range_obj.phases.remove_bad_points(bad_indices)
            range_obj.elevs.remove_bad_points(bad_indices)
            range_obj.remove_bad_alphas(bad_indices)

    @staticmethod
    def filter_inf_lags(range_list): # TODO: Docstring

        for range_obj in range_list:
            log_pwrs = range_obj.pwrs.log_pwrs

            inf_indices = [idx for idx, log_pwr in enumerate(log_pwrs) if not np.isfinite(log_pwr)]

            range_obj.pwrs.remove_bad_points(inf_indices)
            range_obj.remove_bad_alphas(inf_indices)    # REVIEW #0 Shouldn't you have remove bad points for phases/elevs as well?

    @staticmethod
    def filter_low_pwr_lags(raw_data, range_list):
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
        noise_mean = raw_data['noise.mean']  # REVIEW #42 - noise_mean is not used

        # Division by zero error
        if nave <= 0:  # REVIEW #39 - more pythonic to raise exception here, then catch it where you're calling this function
            return

        for range_obj in range_list:
            range_number = range_obj.range_number
            if len(range_obj.pwrs.log_pwrs) == 0:
                continue
            # REVIEW #3 Just a comment about what this actually does, or maybe reference the paper used?
            log_sigma_fluc = np.log(FLUCTUATION_CUTOFF_COEFF *
                                    pwr0[range_number] / math.sqrt(2 * nave))

            bad_indices = []
            cutoff_lag = mplgs + 1

            for idx, (log_pwr, alpha_2) in enumerate(zip(range_obj.pwrs.log_pwrs, # REVIEW #7 Can't seem to find where alpha_2 is defined, but looks like it's called 'power correction factor' in 'Derivations for Phase and Velocity Errors'
                                                         range_obj.pwrs.alpha_2)):
                if idx > cutoff_lag:
                    bad_indices.append(idx)
                else:
                    if((1/np.sqrt(alpha_2) <= ALPHA_CUTOFF) and
                       ((log_pwr < log_sigma_fluc) or np.isclose(log_pwr, log_sigma_fluc))):  # REVIEW refactor? this looks very complex and messy.
                            cutoff_lag = idx
                            bad_indices.append(idx)

            range_obj.pwrs.remove_bad_points(bad_indices) # REVIEW #0 Shouldn't you have remove bad points for phases/elevs and remove_bad_alphas as well?

    @staticmethod
    def filter_bad_acfs(raw_data, ranges, noise_pwr):
        """Removes bad ACFs entirely from analysis

        Removes ACFs which are deemed to be pure noise, or ACFs with too
        few power lags left

        :param raw_data: a dictionary of raw data parameters
        :param range_list: A list of Range objects with data points
        :param noise_pwr: minimum power for which an ACF is pure noise

        """
        nave = raw_data['nave']
        pwr0 = raw_data['pwr0']

        # Division by zero error
        if nave <= 0: # REVIEW #39 - more pythonic to raise exception here, then catch it where you're calling this function
            return

        cutoff_pwr = noise_pwr * 2  # REVIEW Is this just a standard cutoff power? Where is this decided?

        bad_indices = []
        for idx, rang in enumerate(ranges):
            range_number = rang.range_number
            pw = float(pwr0[range_number])  # REVIEW #26 Poor names (pw and ln)
            ln = len(rang.pwrs.log_pwrs)
            if (pw <= cutoff_pwr) or (ln < MIN_LAGS):
                bad_indices.append(idx)
                continue
            # check to see if all powers are equal
            pwr_val = rang.pwrs.log_pwrs[0]  # REVIEW #7 why would you check if all the powers are equal? seems to do nothing, why not just check the len(rang.pwrs.log_pwrs) ?
            for pwr in rang.pwrs.log_pwrs:
                if pwr != pwr_val:
                    break
            else:
                bad_indices.append(idx)

        for i in sorted(bad_indices, reverse=True):
            del ranges[i]

    @staticmethod
    def filter_bad_fits(ranges):  # TODO: Docstring
        bad_indices = []
        for idx, range_obj in enumerate(ranges):
            if range_obj.phase_fit.b == 0.0 or \
               range_obj.linear_pwr_fit.b == 0.0 or \
               range_obj.quadratic_pwr_fit.b == 0.0:
                    bad_indices.append(idx)

        for i in sorted(bad_indices, reverse=True):
            del ranges[i]
