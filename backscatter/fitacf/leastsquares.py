import numpy as np
import math
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

        sigma_2 = sigmas**2

        if fit_type == 'linear':
            lst_sqrs.S = np.sum(np.reciprocal(sigma_2))
            lst_sqrs.S_x = np.sum(x_arr/sigma_2)
            lst_sqrs.S_y = np.sum(y_arr/sigma_2)
            lst_sqrs.S_xx = np.sum((x_arr * x_arr)/sigma_2)
            lst_sqrs.S_xy = np.sum((x_arr * y_arr)/sigma_2)
        elif fit_type == 'quadratic':
            x_2 = x_arr**2
            lst_sqrs.S = np.sum(np.reciprocal(sigma_2))
            lst_sqrs.S_x = np.sum(x_2/sigma_2)
            lst_sqrs.S_y = np.sum(y_arr/sigma_2)
            lst_sqrs.S_xx = np.sum((x_2**2)/sigma_2)
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