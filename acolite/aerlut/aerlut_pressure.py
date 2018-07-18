## def aerlut_pressure
## interpolates aerlut to given pressure
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2017-07-18
## modifications:
##                2018-07-18 (QV) changed acolite import name
def aerlut_pressure(lut, lutdir, pressure, sensor, rsr_file, lut_data_dict=None):
            from acolite.aerlut import get_sensor_lut
        
            lut_split = lut.split('-')
            
            ## get bounding LUT pressures
            lut_pressures = [500,1013,1100] ## typical pressure range for about -500 to +5000m elevation
            lut_sims = ['{}mb'.format(str(lp).zfill(4)) for lp in lut_pressures]
            lut_idx, lut_pressure = (min(enumerate(lut_pressures), key=lambda x: abs(x[1]-pressure)))
            
            lut0_off = 0 if pressure > lut_pressure else -1
            lut1_off = 1 if pressure > lut_pressure else 0

            ## lower bounding LUT
            lut0 = lut_split[0:-1]
            lut0.append(lut_sims[lut_idx+lut0_off])
            lut0 = '-'.join(lut0)

            ## higher bounding LUT
            lut1 = lut_split[0:-1]
            lut1.append(lut_sims[lut_idx+lut1_off])
            lut1 = '-'.join(lut1)
            
            if lut_data_dict is None:
                 ## read sensor LUTs
                 lut0_sensor, meta0_sensor = get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut0, override=0)
                 lut1_sensor, meta1_sensor = get_sensor_lut(sensor, rsr_file, lutdir=lutdir, lutid=lut1, override=0)
            else:
                 lut0_sensor, meta0_sensor = lut_data_dict[lut0]["lut"], lut_data_dict[lut0]["meta"]
                 lut1_sensor, meta1_sensor = lut_data_dict[lut1]["lut"], lut_data_dict[lut1]["meta"]

            ## linear weighting to current pressure
            p0=float(meta0_sensor['press'])
            p1=float(meta1_sensor['press'])
            pc = (pressure - p0) / (p1-p0)
            lut_sensor = {}
            for k in lut0_sensor.keys():
                lut_sensor[k] = (1-pc)*lut0_sensor[k] + pc*lut1_sensor[k]
            
            ## update metadata
            meta_sensor = meta0_sensor.copy()
            meta_sensor['base'] = '-'.join(lut_split[0:-1]+['{}mb'.format(str(pressure).zfill(4)),'interp'])
            meta_sensor['press'] = pressure
            return(lut_sensor, meta_sensor)
