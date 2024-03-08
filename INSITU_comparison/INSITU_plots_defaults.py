xlabel_default = 'HYPSTAR®'
ylabel_default = 'AERONET-OC'


def get_str_stat_list(valid_stats, stat_list, wl):
    str0 = ''
    for stat in stat_list:
        if len(str0) > 0:
            str0 = f'{str0}\n'
        if stat == 'N':
            val = valid_stats[stat]
            str0 = f'{str0}N={val}'
        if stat == 'NMATCH-UPS':
            val = valid_stats['N']
            # self.dataset_w = Dataset(self.path_nc, 'r')
            # nwl = np.unique(np.array(self.dataset_w.variables['mu_wavelength'])).shape[0]
            # self.dataset_w.close()
            # val = val / nwl
            str0 = f'{str0}N={val:.0f}'
        if stat == 'r2':
            val = valid_stats['DETER(r2)']
            str0 = f'{str0}R\u00b2={val:.2f}'
        if stat == 'RMSD' or stat == 'BIAS':
            val = valid_stats[stat]
            if stat == 'BIAS':
                stat = stat.lower()
            # str0 = f'{str0}{stat}={val:.1e}'
            str0 = f'{str0}{stat}={val:.2f}'
        if stat == 'RPD' or stat == 'APD':
            val = valid_stats[stat]
            str0 = f'{str0}{stat}={val:.0f}%'
        if stat == 'WL':
            wls = get_wl_str_from_wl(wl)
            str0 = f'{str0}{wls} nm'
            # str0 = f'{str0}{wl:.2f} nm'

    return str0


def get_wl_str_from_wl(wl_value):
    wl_sat_value_str = f'{wl_value:.2f}'
    if wl_sat_value_str.endswith('.00'):
        return wl_sat_value_str[:-3]
    else:
        return wl_sat_value_str
