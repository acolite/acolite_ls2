## wrapper for acolite CLI/GUI modes
## QV 2018

## import sys to parse arguments
import sys

## fix matplotlib backend (mainly for mac?)
if sys.platform in ['darwin', 'win32']:
    import matplotlib
    matplotlib.use("Agg")

## import acolite source
import acolite as ac

## run command line if --cli provided, otherwise use gui
if '--cli' in sys.argv:
    ac.acolite.acolite_cli(sys.argv)
else:
    ret = ac.acolite.acolite_gui(sys.argv)
