import math

color_dict = dict({ \
    '400.00': 'LightBlue', \
    '412.50': 'DeepSkyBlue', \
    '442.50': 'DodgerBlue', \
    '490.00': 'Blue', \
    '510.00': 'ForestGreen', \
    '560.00': 'Green', \
    '620.00': 'LightCoral', \
    '665.00': 'Red', \
    '673.75': 'Crimson', \
    '681.25': 'FireBrick', \
    '708.75': 'Silver', \
    '753.75': 'Gray', \
    '778.75': 'DimGray', \
    '865.00': 'SlateGray', \
    '885.00': 'DarkSlateGray', \
    '1020.50': 'Pink'})

xlabel_default = r'In situ R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'
ylabel_default = r'Satellite R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'
ylabel_rrs = r'R$_r$$_s$ (sr$^-$$^1$)'
ylabel_rrs_scaled = r'R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'
# xlabel_default = r'Panthyr R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'
# ylabel_default = r'Hypstar R$_r$$_s$ (10$^-$$^3$ sr$^-$$^1$)'

label_insitu_default = 'In situ Rrs'
xlabel_wl_default = 'Wavelength (nm)'

colors_default = ['Blue', 'Red', 'Green', 'Yellow', 'Cyan','Orange','Pink']


def get_color_ref(wlvalue):
    dif_ref = 10000
    color_out = None
    for wlp in color_dict:
        wlpv = float(wlp)
        dif = abs(wlpv - wlvalue)
        if dif < dif_ref:
            dif_ref = dif
            color_out = color_dict[wlp]
    return color_out


def get_color_flag(flagvalue):
    index = int(math.log2(flagvalue))
    return colors_default[index]


def get_color_list(n):
    if n <= len(colors_default):
        return colors_default[0:n]
    else:
        return colors_default[0]
