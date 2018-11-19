## def plot_dark_spectrum
## plots ACOLITE retrieved dark spectrum - to be improved
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-05 (moved from acolite_ac)
## modifications: 2018-06-11 (QV) fixed plotting of S2/MSI results
##                2018-09-12 (QV) fixed matplotlib import issue
##                2018-09-26 (QV) added plt.close() for memory management in long loops

def plot_dark_spectrum(metadata, ds_plot, bands, band_names, data_type, waves, ratm_s, rorayl_s, rdark, rdark_sel, dark_spectrum_option, dark_idx, tau550,sel_model_lut_meta):
                import matplotlib
                import matplotlib.pyplot as plt

                fig = plt.figure()
                canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
                ax = fig.add_subplot(111)

                ratm = []
                rorayl = []
                dark = []
                ds_waves = []
            
                sensor = metadata['SATELLITE_SENSOR']
                sat, sen = sensor.split('_')

                for b, band in enumerate(rdark): #_sorted
                    if data_type == 'NetCDF':
                        btag = '{}'.format(band_names[b].lstrip('B'))

                        if metadata['SENSOR'] == 'MSI':
                            btag = band
                            idx = [i for i,ib in enumerate(band_names) if ib == 'B{}'.format(btag)]
                        else:
                            idx = [i for i,ib in enumerate(band_names) if ib == '{}'.format(band)]
                        
                    if data_type == 'Landsat':
                        if band_names[b] in ['8','9','10','11']: continue
                        idx = [i for i,ib in enumerate(band_names) if ib == '{}'.format(band)]

                    if data_type == 'Sentinel':
                        idx = [i for i,ib in enumerate(band_names) if ib == 'B{}'.format(band)]
                    
                    if len(idx) == 1: 
                        btag = '{}'.format(band_names[idx[0]].lstrip('B'))
                        if btag not in ratm_s: continue
                        ds_waves.append(float(waves[idx[0]]))
                        ratm.append(ratm_s[btag])
                        rorayl.append(rorayl_s[btag])
                        if dark_spectrum_option=='dark_list':
                            dark.append(rdark_sel[btag])
                        else:
                            dark.append(rdark[btag])
                
                if data_type == 'NetCDF':
                    band_title = dark_idx #list(rdark.keys())[dark_idx]
                if data_type == 'Landsat':
                    band_title = dark_idx #bands_sorted[dark_idx]
                if data_type == 'Sentinel':
                    band_title = dark_idx #'B{}'.format(bands_sorted[dark_idx])

                ax.plot(ds_waves, dark, 'o--', color='black', label=r'$\rho_{dark}$')
                ax.plot(ds_waves, rorayl, 'o--', color='blue', label=r'$\rho_{Rayleigh}$')
                ax.plot(ds_waves, ratm, 'o--', color='red', label=r'$\rho_{path}$')
                ax.set_xlim(400,2300)
                ax.set_ylim(0,0.14)
                ax.set_xlabel('Wavelength (nm)')
                ax.set_ylabel(r'$\rho$')

                ax.set_title('{}/{} {}\n{}'.format(sat, sen, metadata['TIME'].strftime('%Y-%m-%d (%H:%M UTC)'),
                                               r'$\theta_s$='+ '{:.1f} '.format(metadata['THS'])+ r'$\tau_{a}550$'+'={:.3f} (mod{}, {})'.format(tau550,sel_model_lut_meta['aermod'][0], band_title)))
                ax.legend(loc='upper right')
                canvas.print_figure(ds_plot, dpi=150)
                plt.close()
                canvas = None
                fig = None
                ax = None
                plt= None
