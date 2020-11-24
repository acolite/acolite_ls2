## def acolite_olci_cli
## simple wrapper for olci processing
## written by Quinten Vanhellemont,
## 2020-11-24
## modifications:

def acolite_olci_cli(*args):

    ## try to import acolite code
    ## find location of current script
    import os, sys
    path = os.path.dirname(os.path.abspath(__file__))
    path_split = (path.split(os.path.sep))

    ## run through the path and try importing acolite code
    imported = False
    for i in range(len(path_split)):
        if path_split[i] == 'acolite':
            ap = os.path.sep.join(path_split[0:i+1])
            if not imported:
                try:
                    sys.path.append(ap)
                    import acolite as ac
                    import acolite.sentinel3
                    imported = True
                except:
                    sys.path.pop()
                    continue
    if not imported:
        print('Could not import acolite source')
        exit()

    ###
    import argparse
    parser = argparse.ArgumentParser(description='Simplified ACOLITE OLCI CLI')
    parser.add_argument('--image', help='Comma separated list of input images')
    parser.add_argument('--output', help='Output directory')
    parser.add_argument('--limit', help='Limit S,W,N,E to crop imagery')

    parser.add_argument('--smile_correction', help='Smile correction (default=True)', default=True)
    parser.add_argument('--use_supplied_pressure', help='Use supplied pressure data (default=True)', default=True)
    parser.add_argument('--use_supplied_altitude', help='Use supplied altitude data (default=True)', default=True)

    parser.add_argument('--tiled_processing', help='Tiled AOT retrieval (default=True)', default=True)
    parser.add_argument('--tile_dimensions', help='AOT tile dimensions in metres', default=[36000,36000])
    parser.add_argument('--write_resolved_atmosphere', help='Write resolved atmosphere retrievals (default=False)', default=False)

    parser.add_argument('--dark_spectrum_option', help='Option to retrieve dark spectrum (default=darkest)', default='darkest')
    parser.add_argument('--verbosity', help='Output verbosity (default=10)', default=10)
    parser.add_argument('--map_rgb', help='Output RGB maps (default=True)', default=True)

    parser.add_argument('--rhod_tgas_cutoff', help='Gas transmittance cut off for dark spectrum (default=0.85)', default=0.85)
    parser.add_argument('--rhod_min_wave', help='Minimum wavelength for dark spectrum (default=440)', default=440.)

    args, unknown = parser.parse_known_args()


    ## for evaluation of bool strings
    import distutils.core

    if args.image is None:
        print('No images given.')
        exit()
    else:
        args.image = args.image.split(',')

    if args.output is None:
        print('No output directory given.')
        exit()

    ## fix limit to 4 element float list if given
    if args.limit is not None:
        limit = [float(s) for s in args.limit.split(',')]
        if len(limit)==4:
            args.limit = limit
        else:
            args.limit = None

    ## check floats
    for k in ['rhod_tgas_cutoff', 'rhod_min_wave']:
        if type(args.__dict__[k]) == str:
            args.__dict__[k] = float(args.__dict__[k])

    ## check bools
    for k in ['smile_correction', 'use_supplied_pressure', 'use_supplied_altitude',
              'tiled_processing', 'write_resolved_atmosphere', 'map_rgb']:
        if type(args.__dict__[k]) == str:
            args.__dict__[k] = bool(distutils.util.strtobool(args.__dict__[k]))

    ## check tile dimensions
    if type(args.tile_dimensions) is str:
        dims = [float(s) for s in args.tile_dimensions.split(',')]
        if len(dims)==2:
            args.tile_dimensions = dims
        else:
            args.limit = [36000,36000]


    print('ACOLITE OLCI CLI settings')
    print(args.__dict__)

    ## set mpl device to agg if plotting rgbs
    if args.map_rgb:
        import matplotlib
        matplotlib.use("Agg")

    ## run processing
    for bundle in args.image:
        ac.sentinel3.acolite_olci(bundle, args.output, limit=args.limit, verbosity=int(args.verbosity),
                                smile_correction=args.smile_correction,
                                use_supplied_pressure=args.use_supplied_pressure,
                                use_supplied_altitude=args.use_supplied_altitude,
                                dark_spectrum_option=args.dark_spectrum_option,
                                rhod_tgas_cutoff=args.rhod_tgas_cutoff,
                                rhod_min_wave=args.rhod_min_wave,
                                tiled_processing=args.tiled_processing,
                                write_resolved_atmosphere = args.write_resolved_atmosphere,
                                tile_dimensions = args.tile_dimensions,
                                map_rgb = args.map_rgb)

if __name__ == '__main__':
    acolite_olci_cli()
