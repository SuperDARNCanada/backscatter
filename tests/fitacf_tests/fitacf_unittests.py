import backscatter.fitacf as fit
import unittest
import backscatter.dmap as dm

class FakeRangeClassForTesting(fit.Range):
    def __init__(self, range_number):
        self.range_number = range_number

raw_data = dm.parse_dmap_format_from_file('20131004.0401.00.rkn.rawacf')
test_range_obj = FakeRangeClassForTesting(0)
lags = fit.create_lag_list(raw_data[0])

class TestFitacf(unittest.TestCase):

	def test_lag_list(self):
		print("\n\ntesting creating lag list")
		lags = fit.create_lag_list(raw_data[0])
		print(lags)

	def test_find_CRI(self):
		print("\n\nStarting test_find_CRI")
		global test_range_obj
		test_range_obj.CRI = test_range_obj.find_cri(raw_data[0])
		print(test_range_obj.CRI)

	def test_find_alphas(self):
		print("\n\nStarting test_find_alphas")
		global test_range_obj
		#lags = fit.create_lag_list(raw_data[0])
		test_range_obj.alpha_2 = test_range_obj.find_alphas(raw_data[0],lags)

		print(test_range_obj.alpha_2)

	def test_create_powers(self):
		print("\n\nStarting test_create_powers")
		#lags = fit.create_lag_list(raw_data[0])
		pwr_data_points = fit.PowerDataPoints(raw_data[0],lags,test_range_obj)
		print("power points")
		print(pwr_data_points.log_pwrs)
		print("power sigmas")
		print(pwr_data_points.sigmas)
		print("power t values")
		print(pwr_data_points.t)

	def test_create_acf_phases(self):
		print("\n\nStarting test_create_acf_phases")
		global test_range_obj
		#lags = fit.create_lag_list(raw_data[0])
		test_range_obj.CRI = test_range_obj.find_cri(raw_data[0])
		test_range_obj.alpha_2 = test_range_obj.find_alphas(raw_data[0],lags)
		phase_data_points = fit.PhaseDataPoints(raw_data[0],'acfd',lags,test_range_obj)
		print("phase points")
		print(phase_data_points.phases)
		print("phase sigmas")
		print(phase_data_points.sigmas)
		print("phase t values")
		print(phase_data_points.t)

	def test_create_xcf_phases(self):
		print("\n\nStarting test_create_xcf_phases")
		global test_range_obj
		test_range_obj.CRI = test_range_obj.find_cri(raw_data[0])
		test_range_obj.alpha_2 = test_range_obj.find_alphas(raw_data[0],lags)
		phase_data_points = fit.PhaseDataPoints(raw_data[0],'xcfd',lags,test_range_obj)
		print("phase points")
		print(phase_data_points.phases)
		print("phase sigmas")
		print(phase_data_points.sigmas)
		print("phase t values")
		print(phase_data_points.t)

	def test_filter_tx_overlap(self):
		print("\n\nStarting test_filter_tx_overlap")
		test_range_obj_2 = fit.Range(20,raw_data[0],lags)
		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,[test_range_obj_2])
		print("power points")
		print(test_range_obj_2.pwrs.log_pwrs)

	def test_filter_low_pwr_lags(self):
		print("\n\nStarting test_filtering_low_pwr_lags")
		test_range_obj_list = []
		for i in range(0,raw_data[0]['nrang']):
			test_range_obj_list.append(fit.Range(i,raw_data[0],lags))

		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,test_range_obj_list)


		fit.Filtering.filter_low_pwr_lags(raw_data[0],test_range_obj_list)

		counts = []
		print("counts:")
		for tr in test_range_obj_list:
			counts.append(len(tr.pwrs.log_pwrs))

		print(counts)

		print("range {0} powers".format(2))
		print(test_range_obj_list[2].pwrs.log_pwrs)

	def test_filter_acfs(self):
		print("\n\nStarting test_filtering_bad_acfs")
		test_range_obj_list = []
		# for data in raw_data:
		# 	for i in range(0,data['nrang']):
		# 		if data['pwr0'][i] != 0:
		# 			test_range_obj_list.append(fit.Range(i,data,lags))

		for i in range(0,raw_data[0]['nrang']):
			if raw_data[0]['pwr0'][i] != 0:
				test_range_obj_list.append(fit.Range(i,raw_data[0],lags))

		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,test_range_obj_list)


		fit.Filtering.filter_low_pwr_lags(raw_data[0],test_range_obj_list)

		noise_pwr = fit.NoisePower.acf_cutoff_pwr(raw_data[0])
		fit.Filtering.filter_bad_acfs(raw_data[0],test_range_obj_list,noise_pwr)

		ranges_left = []
		print("ranges left:")
		for tr in test_range_obj_list:
			ranges_left.append(tr.range_number)

		print(ranges_left)

	def test_power_fitting(self):
		print("\n\nStarting test_power_fitting")
		test_range_obj_list = []
		# for data in raw_data:
		# 	for i in range(0,data['nrang']):
		# 		if data['pwr0'][i] != 0:
		# 			test_range_obj_list.append(fit.Range(i,data,lags))

		for i in range(0,raw_data[0]['nrang']):
			if raw_data[0]['pwr0'][i] != 0:
				test_range_obj_list.append(fit.Range(i,raw_data[0],lags))

		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,test_range_obj_list)


		fit.Filtering.filter_low_pwr_lags(raw_data[0],test_range_obj_list)

		noise_pwr = fit.NoisePower.acf_cutoff_pwr(raw_data[0])
		fit.Filtering.filter_bad_acfs(raw_data[0],test_range_obj_list,noise_pwr)

		fit.ACFFitting.acf_pwr_fitting(test_range_obj_list)

		print(vars(test_range_obj_list[0].linear_pwr_fit))

	def test_calculate_sigmas_for_phase(self):
		pass

	def test_acf_phase_unwrap(self):
		pass

	def test_acf_phase_fit(self):
		print("\n\nStarting test_acf_phase_fitting")
		test_range_obj_list = []
		# for data in raw_data:
		# 	for i in range(0,data['nrang']):
		# 		if data['pwr0'][i] != 0:
		# 			test_range_obj_list.append(fit.Range(i,data,lags))

		for i in range(0,raw_data[0]['nrang']):
			if raw_data[0]['pwr0'][i] != 0:
				test_range_obj_list.append(fit.Range(i,raw_data[0],lags))

		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,test_range_obj_list)


		fit.Filtering.filter_low_pwr_lags(raw_data[0],test_range_obj_list)

		noise_pwr = fit.NoisePower.acf_cutoff_pwr(raw_data[0])
		fit.Filtering.filter_bad_acfs(raw_data[0],test_range_obj_list,noise_pwr)
		fit.ACFFitting.acf_pwr_fitting(test_range_obj_list)

		fit.ACFFitting.calculate_phase_and_elev_sigmas(test_range_obj_list,raw_data[0])
		fit.ACFFitting.acf_phase_unwrap(test_range_obj_list)
		fit.ACFFitting.acf_phase_fitting(test_range_obj_list)

		print(vars(test_range_obj_list[0].phase_fit))		

	def test_xcf_phase_fit(self):
		print("\n\nStarting test_xcf_phase_fitting")
		test_range_obj_list = []
		# for data in raw_data:
		# 	for i in range(0,data['nrang']):
		# 		if data['pwr0'][i] != 0:
		# 			test_range_obj_list.append(fit.Range(i,data,lags))

		for i in range(0,raw_data[0]['nrang']):
			if raw_data[0]['pwr0'][i] != 0:
				test_range_obj_list.append(fit.Range(i,raw_data[0],lags))


		fit.Filtering.filter_tx_overlapped_lags(raw_data[0],lags,test_range_obj_list)


		fit.Filtering.filter_low_pwr_lags(raw_data[0],test_range_obj_list)

		noise_pwr = fit.NoisePower.acf_cutoff_pwr(raw_data[0])
		fit.Filtering.filter_bad_acfs(raw_data[0],test_range_obj_list,noise_pwr)
		fit.ACFFitting.acf_pwr_fitting(test_range_obj_list)

		fit.ACFFitting.calculate_phase_and_elev_sigmas(test_range_obj_list,raw_data[0])
		fit.ACFFitting.acf_phase_unwrap(test_range_obj_list)
		fit.ACFFitting.acf_phase_fitting(test_range_obj_list)

		fit.ACFFitting.xcf_phase_unwrap(test_range_obj_list)
		fit.ACFFitting.xcf_phase_fitting(test_range_obj_list)

		print(test_range_obj_list[0].range_number,test_range_obj_list[0].elevs.sigmas)
		print(vars(test_range_obj_list[0].elev_fit))	





if __name__ == '__main__':
	unittest.main()