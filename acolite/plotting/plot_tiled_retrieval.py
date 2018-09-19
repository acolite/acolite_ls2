## def plot_tiled_retrieval
## plots ACOLITE tiled aot retrievals
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-05
## modifications: 2018-09-19 (QV) changed matplotlib imports

def plot_tiled_retrieval(tiles, tile_data, tile_output, metadata, rdark, odir,
                        map_model_per_tile=True, map_rdark_blue_swir=True, map_aot_blue_swir=True):
    import matplotlib.cm
    import matplotlib.figure
    import matplotlib.backends.backend_agg

    from numpy.ma import masked_where
    from numpy import isnan
    
    cmap=matplotlib.cm.viridis

    sensor = metadata['SATELLITE_SENSOR']
    sat, sen = sensor.split('_')

    ### map aot, band and model per tile
    if map_model_per_tile:
        band_labels=list(rdark.keys())
        band_ticks=[b+1 for b,bl in enumerate(band_labels)]
 
        fig = matplotlib.figure.Figure(figsize=(16,5))
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)

        ax1 = fig.add_subplot(1,3,1)
        ax2 = fig.add_subplot(1,3,2)
        ax3 = fig.add_subplot(1,3,3)

        im=ax1.imshow(tile_output['tau550'], vmin=0, vmax=0.5)
        fig.colorbar(im, ax=ax1,orientation='horizontal', label=r'$\tau_a$ 550')

        im=ax2.imshow(tile_output['band'], vmin=min(band_ticks), vmax=max(band_ticks))
        cbar = fig.colorbar(im, ax=ax2,orientation='horizontal', ticks=band_ticks, label='band')
        cbar.ax.set_xticklabels(band_labels)

        im=ax3.imshow(tile_output['model'], vmin=1, vmax=3)
        cbar = fig.colorbar(im, ax=ax3,orientation='horizontal', ticks=[1,2,3], label='model')
        cbar.ax.set_xticklabels(['C','M','U'])

        #for a in ax: a.axis('off')
        ax1.axis('off')
        ax2.axis('off')
        ax3.axis('off')

        plfile = '{}/{}.png'.format(odir, '{}_{}_AOT_tiled_{}x{}'.format(metadata['SATELLITE_SENSOR'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tiles[0],tiles[1]))
        canvas.print_figure(plfile, dpi=300)
        
    ### map rdark for two bands (blue and swir)
    if map_rdark_blue_swir:
        band_labels=list(rdark.keys())
        band_ticks=[b+1 for b,bl in enumerate(band_labels)]

        fig = matplotlib.figure.Figure(figsize=(12,6))
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)

        k1=list(tile_data.keys())[0]
        k2=list(tile_data.keys())[-1 if sen == 'MSI' else -4]

        im=ax1.imshow(tile_data[k1]['rdark'], vmin=0, vmax=0.25, cmap=cmap)
        cbar= fig.colorbar(im, ax=ax1,orientation='horizontal', label=r'$\rho_{dark}$'+' band {}'.format(k1))

        im=ax2.imshow(tile_data[k2]['rdark'], vmin=0, vmax=0.05, cmap=cmap)
        cbar= fig.colorbar(im, ax=ax2,orientation='horizontal', label=r'$\rho_{dark}$'+' band {}'.format(k2))

        #for a in ax: a.axis('off')
        ax1.axis('off')
        ax2.axis('off')

        plfile = '{}/{}.png'.format(odir, '{}_{}_rdark_tiled_{}x{}'.format(metadata['SATELLITE_SENSOR'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tiles[0],tiles[1]))
        canvas.print_figure(plfile, dpi=300)


    ### map aot retrieved from two bands (blue and swir)
    if map_aot_blue_swir:
        band_labels=list(rdark.keys())
        band_ticks=[b+1 for b,bl in enumerate(band_labels)]

        fig = matplotlib.figure.Figure(figsize=(12,6))
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        ax1 = fig.add_subplot(1,2,1)
        ax2 = fig.add_subplot(1,2,2)

        k1=list(tile_data.keys())[0]
        k2=list(tile_data.keys())[-1 if sen == 'MSI' else -4]

        im=ax1.imshow(tile_data[k1]['tau550'], vmin=0, vmax=0.5, cmap=cmap)
        cbar= fig.colorbar(im, ax=ax1,orientation='horizontal', label=r'$\tau_{a} 550$'+' (band {})'.format(k1))

        im=ax2.imshow(tile_data[k2]['tau550'], vmin=0, vmax=0.5, cmap=cmap)
        cbar= fig.colorbar(im, ax=ax2,orientation='horizontal', label=r'$\tau_{a} 550$'+' (band {})'.format(k2))

        #for a in ax: a.axis('off')
        ax1.axis('off')
        ax2.axis('off')

        plfile = '{}/{}.png'.format(odir, '{}_{}_tauband_tiled_{}x{}'.format(metadata['SATELLITE_SENSOR'],metadata['TIME'].strftime('%Y_%m_%d_%H_%M_%S'), tiles[0],tiles[1]))
        canvas.print_figure(plfile, dpi=300)
