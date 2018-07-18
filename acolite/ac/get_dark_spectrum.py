## def get_dark_spectrum
## gets dark spectrum according to percentiles
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-12-12
## modifications: 2016-12-13 (QV) added rorayl option to select pixels as least as bright as the Rayleigh reflectance
##                2017-01-23 (QV) added absolute pixel option
##                2017-01-19 (QV) added finite check for absolute pixel
##                2017-11-21 (QV) added dark_list alias
##                2018-02-08 (QV) fixed max_len for dark_list option (doesn't work band per band)
##                2018-03-07 (QV) changed 'is' to '=='

def get_dark_spectrum(data, option='percentile', rorayl=None, percentile=0.1, 
                     pixel_idx=10, perc_idx= 1, percentiles = [0,0.1,1,5,10,25,50,75,90,95,99,99.9,100],
                     pixel_range_min=0, pixel_range_max=1000):
            
            if (option == 'percentile'):
                from numpy import nanpercentile
                ## get percentiles and set up rtoa dark spectrum
                all_perc = {}
                dark_spectrum = {}
                for band in data.keys():
                    all_perc[band] = nanpercentile(data[band], percentiles)
                    dark_spectrum[band] = all_perc[band][perc_idx]

            if (option == 'minRayleigh') & (rorayl is not None):
                from numpy import nanpercentile
                all_perc = {}
                dark_spectrum = {}
                for band in data.keys():
                    all_perc[band] = nanpercentile(data[band], percentiles)
                    sub = data[band][data[band] > rorayl[band]]
                    dark_spectrum[band] = nanpercentile(sub, percentile)

            if (option == 'sortRayleigh') & (rorayl is not None):
                from numpy import sort, nanpercentile
                all_perc = {}
                dark_spectrum = {}
                for band in data.keys():
                    all_perc[band] = nanpercentile(data[band], percentiles)
                    sub = sort(data[band][data[band] > rorayl[band]], axis=None)[0:500]
                    dark_spectrum[band] = nanpercentile(sub, 1.)

            if (option == 'absolute_pixel'):
                from numpy import sort, nanpercentile, isfinite, nan
                all_perc = {}
                dark_spectrum = {}
                for band in data.keys():
                    all_perc[band] = nanpercentile(data[band], percentiles)
                    dsorted = sort(data[band][isfinite(data[band])], axis=None)
                    if len(dsorted) > pixel_idx:
                        dark_spectrum[band] = dsorted[pixel_idx]
                    else:
                        dark_spectrum[band] = nan

            if (option == 'absolute_pixel_list') or (option == 'dark_list'):
                from numpy import sort, nanpercentile, isfinite, nan
                all_perc = {}
                dark_spectrum = {}
                max_len = 1000
                max_len = pixel_range_max - pixel_range_min
                #min_len = max_len * 1.0
                for band in data.keys():
                    all_perc[band] = nanpercentile(data[band], percentiles)
                    dsorted = sort(data[band][isfinite(data[band])], axis=None)
                    dark_spectrum[band] = dsorted[pixel_range_min:pixel_range_max]
                    if max_len > len(dark_spectrum[band]): max_len = len(dark_spectrum[band])
                    #if min_len > len(dark_spectrum[band]): min_len = len(dark_spectrum[band])
                for band in data.keys():
                    #print(band, len(dark_spectrum[band]), max_len)
                    if len(dark_spectrum[band]) > max_len: 
                        dark_spectrum[band] = dark_spectrum[band][0:max_len]

            return dark_spectrum, all_perc
