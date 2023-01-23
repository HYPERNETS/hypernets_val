
class SAT_EXTRACTS_LIST:
    def __init__(self, boptions, verbose):
        self.boptions = boptions

    def get_list_as_dict(self):
        start_date = self.boptions.stat_date
        end_date = self.boptions.end_date
        sat_list = {}
        
        # self.param_sat = {
        #     'satellite': sat_satellite.upper(),
        #     'sensor': sat_sensor.upper(),
        #     'platform': sat_platform.upper(),
        #     'resolution': sat_res.upper(),
        #     'ac': atm_corr.upper(),
        #     'prefix': prefix
        # }