import math
import numpy as np

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