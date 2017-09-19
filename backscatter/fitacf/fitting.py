import leastsquares as ls
import numpy as np
import math


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
        lst_sqrs_fit = ls.LeastSquaresFitting(1,1)
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
        lst_sqrs_fit = ls.LeastSquaresFitting(1,1)
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
        lst_sqrs_fit = ls.LeastSquaresFitting(1,1)
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

            phase_inverse_alpha_2 = 1/phases.alpha_2
            elev_inverse_alpha_2 = 1/elevs.alpha_2

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
                #eprint(error_string,phase_inverse_alpha_2,pwr_values)

            if np.isnan(elev_sigmas).any() or np.isinf(elev_sigmas).any():
                error_string = "Elevation sigmas bad at range {0} -- elev_inverse_alphas,pwr_values"
                error_string.format(range_obj.range_number)
                #eprint(error_string,elev_inverse_alpha_2,pwr_values)

            phase_sigmas = np.array([math.pi if sigma > math.pi else sigma for sigma in phase_sigmas])
            elev_sigmas = np.array([math.pi if sigma > math.pi else sigma for sigma in elev_sigmas])

            """Since lag 0 phase is included for elevation fit, we set lag 0 sigma the
            same as lag 1 sigma"""
            elev_sigmas[0] = elev_sigmas[1]

            phases.set_sigmas(phase_sigmas)
            elevs.set_sigmas(elev_sigmas)

    @staticmethod
    def acf_phase_unwrap(range_list,raw_data):
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

            slope_numerator, slope_denominator = 0, 0

            #This is to skip the first element in lists
            phase_prev = phase_values[0]
            sigma_prev = sigma_values[0]
            t_prev = t_values[0]
            orig_phase_iterator = np.nditer([phase_values,sigma_values,t_values])
            orig_phase_iterator.iternext()
            while not orig_phase_iterator.finished:
                phase_curr = orig_phase_iterator[0]
                sigma_curr = orig_phase_iterator[1]
                t_curr = orig_phase_iterator[2]

                phase_diff = phase_curr - phase_prev

                sigma_bar = (sigma_curr + sigma_prev)/2

                t_diff = t_curr - t_prev

                if np.fabs(phase_diff) < math.pi:
                    slope_numerator = slope_numerator + phase_diff/(sigma_bar**2)/t_diff
                    slope_denominator = slope_denominator + 1/(sigma_bar**2)

                phase_prev = phase_curr
                sigma_prev = sigma_curr
                t_prev = t_curr

                orig_phase_iterator.iternext()


            piecewise_slope_estimate = slope_numerator / slope_denominator
            new_phases, total_2pi_corrections = phase_correction(piecewise_slope_estimate,
                                                                    phase_values,t_values)

            if total_2pi_corrections > 0:

                corr_phase_iterator = np.nditer([new_phases,sigma_values,t_values])
                S_xx, S_xy = 0.0, 0.0
                while not corr_phase_iterator.finished:
                    phase = corr_phase_iterator[0]
                    sigma = corr_phase_iterator[1]
                    t = corr_phase_iterator[2]
                    if sigma > 0.0:
                        S_xy = S_xy + (phase * t)/(sigma**2)
                        S_xx = S_xx + (t**2)/(sigma**2)

                    corr_phase_iterator.iternext()

                corr_slope_estimate = S_xy / S_xx

                corr_phase_iterator.reset()
                corr_slope_err = 0.0;
                while not corr_phase_iterator.finished:
                    phase = corr_phase_iterator[0]
                    sigma = corr_phase_iterator[1]
                    t = corr_phase_iterator[2]
                    if sigma > 0.0:
                        corr_slope_err = corr_slope_err + (corr_slope_estimate * t - phase) ** 2 / sigma ** 2

                    corr_phase_iterator.iternext()

                orig_phase_iterator.reset()
                S_xx, S_xy = 0.0, 0.0
                while not orig_phase_iterator.finished:
                    phase = orig_phase_iterator[0]
                    sigma = orig_phase_iterator[1]
                    t = orig_phase_iterator[2]
                    if sigma > 0.0:
                        S_xy = S_xy + (phase * t)/(sigma**2)
                        S_xx = S_xx + (t**2)/(sigma**2)

                    orig_phase_iterator.iternext()

                orig_slope_est = S_xy / S_xx

                orig_phase_iterator.reset()
                orig_slope_err = 0.0;
                while not orig_phase_iterator.finished:
                    phase = orig_phase_iterator[0]
                    sigma = orig_phase_iterator[1]
                    t = orig_phase_iterator[2]
                    if sigma > 0.0:
                        orig_slope_err = (orig_slope_err + (orig_slope_est * t - phase) ** 2
                                            / sigma ** 2)

                    orig_phase_iterator.iternext()


                if (orig_slope_err > corr_slope_err):
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

            new_phases = phase_correction(range_obj.phase_fit.b,elev_phase_values,t_values)[0]

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

            slope_estimate = np.divide(S_xy , S_xx)
            new_phases = phase_correction(slope_estimate,new_phases,t_values)[0]

            range_obj.elevs.set_phases(new_phases)

    @staticmethod
    def phase_correction(slope_estimate,phase_values,t_values):
        """Adds the estimated number of 2*pi corrections to phase values

        :param slope_estimate: predicted slope from iterative unwrap
        :param phase_value: a particular phase to correct
        :param t_value: a value in time for a particular phase point
        :returns: phase shifted by number 2*pi corrections, total 2pi corrections.

        """


        phase_predicted = slope_estimate * t_values

        #I add a rounding here so that if there is inexact division of pi/2pi then
        #.49999999... gets rounded up first.
        phase_diff = np.around((phase_predicted - phase_values)/(2 * math.pi), decimals=5)
        phase_correction = np.array([round(pd) for pd in phase_diff])

        corrected_phase = phase_values + (phase_correction * 2 * math.pi)
        total_corrections = np.sum(np.abs(phase_correction))
        return corrected_phase, total_corrections