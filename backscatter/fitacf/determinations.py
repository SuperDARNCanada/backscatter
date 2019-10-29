from datetime import datetime
import numpy as np
import math
import sys

from backscatter import config
from backscatter import hdw_info as hdw

V_MAX = float(config.get("fitacf","v_max"))
W_MAX = float(config.get("fitacf","w_max"))
C = 299792458.0


FITACF_REVISION_MAJOR = int(config.get("fitacf","fitacf_revision_major"))
FITACF_REVISON_MINOR = int(config.get("fitacf","fitacf_revision_minor"))

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

        #if the range list is empty at this point, no other fields are to be written.
        if not range_list:
            return new_parameter_dict

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

        number_of_good_data = len(new_parameter_dict['qflg'])
        float_zeroes = np.zeros(number_of_good_data)
        int_zeroes = float_zeroes.astype(int)

        new_parameter_dict['x_qflg'] = int_zeroes
        new_parameter_dict['x_gflg'] = int_zeroes
        new_parameter_dict['x_p_l'] = float_zeroes
        new_parameter_dict['x_p_l_e'] = float_zeroes
        new_parameter_dict['x_p_s'] = float_zeroes
        new_parameter_dict['x_p_s_e'] = float_zeroes
        new_parameter_dict['x_v'] = float_zeroes
        new_parameter_dict['x_v_e'] = float_zeroes
        new_parameter_dict['x_w_l'] = float_zeroes
        new_parameter_dict['x_w_l_e'] = float_zeroes
        new_parameter_dict['x_w_s'] = float_zeroes
        new_parameter_dict['x_w_s_e'] = float_zeroes

        if 'xcfd' not in raw_data:
            new_parameter_dict['phi0'] = float_zeroes
        else:
            new_parameter_dict['phi0'] = self.set_xcf_phi0(range_list,raw_data, hdw_info)

        new_parameter_dict['phi0_e'] = self.set_xcf_phi0_err(range_list)

        if 'xcfd' not in raw_data:
            elv = {'low' : float_zeroes, 'normal' : float_zeroes, 'high' : float_zeroes}
        else:
            elv = self.find_elevation(range_list,raw_data,hdw_info)

        new_parameter_dict['elv_low'] = elv['low']
        new_parameter_dict['elv'] = elv['normal']
        new_parameter_dict['elv_high'] = elv['high']
        new_parameter_dict['x_sd_l'] = float_zeroes
        new_parameter_dict['x_sd_s'] = float_zeroes
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

        noise_dB = 10 * np.log10(noise_pwr);

        p_l_conversion = lambda x: 10 * x/np.log(10) - noise_dB

        p_l = [p_l_conversion(range_obj.linear_pwr_fit.a) for range_obj in range_list]

        return np.array(p_l)


    def set_p_l_err(self,range_list):
        """Computes the power in dB of the linear power fit error at each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted lambda power errors in dB

        """

        p_l_err_conversion = lambda x: 10 * np.sqrt(x) / np.log(10)

        p_l_err = [p_l_err_conversion(range_obj.linear_pwr_fit_err.sigma_2_a) for range_obj in range_list]

        return np.array(p_l_err)


    def set_p_s(self,range_list,noise_pwr):
        """Computes the power in dB of the quadratic power fit at each range

        :param range_list: a list of Range objects after fitting
        :param noise_pwr: minimum power for which an ACF is pure noise
        :returns: an array of fitted sigma powers in dB

        """

        noise_dB = 10 * np.log10(noise_pwr)

        p_s_conversion = lambda x: 10 * x/np.log(10) - noise_dB

        p_s = [p_s_conversion(range_obj.quadratic_pwr_fit.a) for range_obj in range_list]

        return np.array(p_s)


    def set_p_s_err(self,range_list):
        """Computes the power in dB of the quadratic power fit error at each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted sigma power errors in dB

        """

        p_s_err_conversion = lambda x: 10 * np.sqrt(x) / np.log(10)

        p_s_err = [p_s_err_conversion(range_obj.quadratic_pwr_fit_err.sigma_2_a) for range_obj in range_list]

        return np.array(p_s_err)


    def set_v(self,range_list,raw_data,hdw_info):
        """Computes the fitted velocity at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :param hdw_info: a dictionary of the radar hardware information
        :returns: an array of computed velocites

        """

        vel_conversion = C/((4*np.pi)*(raw_data['tfreq'] * 1000.0)) * hdw_info['velsign']

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
        vel_conversion = C/((4*np.pi)*(raw_data['tfreq'] * 1000.0)) * hdw_info['velsign']

        vel_err = [np.sqrt(range_obj.phase_fit.sigma_2_b) * vel_conversion for range_obj in range_list]

        return np.array(vel_err)


    def set_w_l(self,range_list,raw_data):
        """Computes the spectral width from the linear power fit at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral widths

        """

        width_conversion = C/((4*np.pi)*(raw_data['tfreq'] * 1000.0))*2.0

        w_l_calculation = lambda x: np.fabs(x) * width_conversion
        w_l = [w_l_calculation(range_obj.linear_pwr_fit.b) for range_obj in range_list]

        return np.array(w_l)


    def set_w_l_err(self,range_list,raw_data):
        """Computes the spectral width error from the linear power fit errors
        at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral width errors

        """

        width_conversion = C/((4*np.pi)*(raw_data['tfreq'] * 1000.0))*2.0

        w_l_err_calculation = lambda x: np.sqrt(x) * width_conversion
        w_l_err = [w_l_err_calculation(range_obj.linear_pwr_fit_err.sigma_2_b) for range_obj in range_list]

        return np.array(w_l_err)

    def set_w_s(self,range_list,raw_data):
        """Computes the spectral width from the quadratic power fit at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral widths

        """

        w_s_conversion = C/(4*np.pi)/(raw_data['tfreq'] * 1000.0) *4.* np.sqrt(np.log(2))

        w_s_calculation = lambda x: np.sqrt(np.fabs(x)) * w_s_conversion
        w_s = [w_s_calculation(range_obj.quadratic_pwr_fit.b) for range_obj in range_list]

        return np.array(w_s)

    def set_w_s_err(self,range_list,raw_data):
        """Computes the spectral width error from the quadratic power fit errors
        at each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :returns: an array of computed spectral width errors

        """

        w_s_conversion = C/(4*np.pi)/(raw_data['tfreq'] * 1000.0) *4.* np.sqrt(np.log(2))

        w_s_calculation = lambda x,y: np.sqrt(x)/2./np.sqrt(np.fabs(y)) * w_s_conversion

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

        antenna_sep = np.sqrt(x*x + y*y + z*z)

        elev_corr = np.arcsin(z/antenna_sep)

        elevations = {}

        if y > 0.0:
            phi_sign = 1
        else:
            phi_sign = -1
            elev_corr = -elev_corr

        azi_offset = hdw_info['maxbeams']/2 - 0.5
        phi_0 = hdw_info['beamsep'] * ( bmnum- azi_offset) * np.pi/180
        c_phi_0 = np.cos(phi_0)

        wave_num = 2 * np.pi * tfreq * 1000/C;

        cable_offset = -2 * np.pi * tfreq * 1000 * hdw_info['tdiff'] * 1.0e-6

        phase_diff_max = phi_sign * wave_num * antenna_sep * c_phi_0 + cable_offset

        psi_calculation = lambda x: x + 2 * np.pi * np.floor((phase_diff_max-x)/(2*np.pi))
        psi_uncorrected = [psi_calculation(range_obj.elev_fit.a) for range_obj in range_list]

        if(phi_sign < 0):
            psi_uncorrected = [psi_u + 2 * np.pi for psi_u in psi_uncorrected]

        psi = [psi_u - cable_offset for psi_u in psi_uncorrected]

        psi_kd = [p/(wave_num * antenna_sep) for p in psi]
        theta = [c_phi_0**2 - pkd**2 for pkd in psi_kd]

        elev_calculation = lambda x: -elev_corr if (x < 0.0 or np.fabs(x) > 1.0) else np.arcsin(np.sqrt(x))
        elevation = [elev_calculation(t) for t in theta]

        elevations['high'] = [180/np.pi * (elev + elev_corr) for elev in elevation]

        #Elevation errors
        psi_k2d2 = [p/(wave_num**2 * antenna_sep**2) for p in psi]
        df_by_dy = [pkd/np.sqrt(t * (1 - t)) for pkd,t in zip(psi_k2d2,theta)]

        elev_low_calculation = lambda x,y: 180/np.pi * np.sqrt(x) * np.fabs(y)
        errors = [range_obj.elev_fit.sigma_2_a for range_obj in range_list]
        elevations['low'] = [elev_low_calculation(err,dfdy) for err,dfdy in zip(errors,df_by_dy)]

        #Experiment to compare fitted and measured elevation
        xcfd = raw_data['xcfd']
        real = [xcfd[range_obj.range_idx][0][0] for range_obj in range_list]
        imag = [xcfd[range_obj.range_idx][0][1] for range_obj in range_list]
        xcf0_p = [np.arctan2(i,r) for i,r in zip(imag,real)]

        psi_uu_calculation = lambda x: x + 2 * np.pi * np.floor((phase_diff_max-x)/(2*np.pi))
        psi_uncorrected_unfitted = [psi_uu_calculation(x) for x in xcf0_p]

        psi_uu_calculation = lambda x: x + (2 * np.pi) if phi_sign < 0 else x
        psi_uncorrected_unfitted = [psi_uu_calculation(p_uu) for p_uu in psi_uncorrected_unfitted]

        psi = [p_uu - cable_offset for p_uu in psi_uncorrected_unfitted]

        psi_kd = [p/(wave_num * antenna_sep) for p in psi]
        theta = [c_phi_0**2 - pkd**2 for pkd in psi_kd]

        elev_calculation = lambda x: -elev_corr if (x < 0.0 or np.fabs(x) > 1.0) else np.arcsin(np.sqrt(x))
        elevation = [elev_calculation(t) for t in theta]

        elevations['normal'] = [180/np.pi * (elev + elev_corr) for elev in elevation]

        return elevations


    def set_xcf_phi0(self, range_list, raw_data, hdw_info):
        """Sets the unfitted offset of the XCF phase for each range

        :param range_list: a list of Range objects after fitting
        :param raw_data: a dictionary of raw data parameters
        :param hdw_info: a dictionary of the radar hardware information
        :returns: an array of unfitted XCF phases offsets

        """
        #phi0 = [range_obj.elev_fit.a for range_obj in range_list]
        xcfd = [raw_data['xcfd'][range_obj.range_idx] for range_obj in range_list]
        phi0 = [np.arctan2(xcf[0][1],xcf[0][0]) for xcf in xcfd]
        phi0 = [p * hdw_info['phasesign'] for p in phi0]

        return np.array(phi0)


    def set_xcf_phi0_err(self,range_list):
        """Sets the fitted offset error of the XCF phase for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of fitted XCF phases offset errors

        """
        phi0_err = [np.sqrt(range_obj.elev_fit.sigma_2_a) for range_obj in range_list]

        return np.array(phi0_err)

    def set_xcf_sdev_phi(self,range_list):
        """Sets the chi_2 value for the XCF phase fit for each range

        :param range_list: a list of Range objects after fitting
        :returns: an array of chi_2 values

        """
        sdev_phi = [range_obj.elev_fit.chi_2 for range_obj in range_list]

        return np.array(sdev_phi)
