## def acolite_l2w_qaa
## runs qaa on dataset read from l2r ACOLITE output file
##
## Written by Quinten Vanhellemont 2018-03-08
## Last modifications:
##                2018-07-18 (QV) changed acolite import name

def acolite_l2w_qaa(data_read, par_attributes, ths = 0.0, satellite = None, compute_zeu_lee = False, mask = None, 
                    k_atag = 'a490_QAAsw', k_bbptag = 'bbp490_QAAsw'):
    from numpy import pi, where, sqrt, power, nan, log10, exp, log, sin, cos, roots, isfinite, real
    from acolite.shared import coef_qaa
    import time

    def qaa_kd(a, bbp, coef):
        from numpy import exp
        v=coef['m'][1]*(1.-coef['m'][2]*exp(-coef['m'][3]*a))
        kd=coef['m'][0]*a+v*bbp
        return(kd)

    coef = coef_qaa()
    qaa_wave = [443, 490, 560, 665]
    sen_wave = par_attributes['waves']
    qaa_data = {}

    if satellite is not None:
        sat = satellite[0:2]
    else:
        if coef['spectral_shift'] == '1':
            print('QAA: No satellite given, not performing spectral shift.')
            coef['spectral_shift'] = '0'

    print('Started QAA computations')
    t0 = time.time()
    ## get data and convert to subsurface rrs
    for iw,wave in enumerate(sen_wave):
        rrs = data_read['{}_{}'.format(par_attributes['dataset'],wave)]/pi

        ## do band shifting with coefficients provided by Dimitry
        if coef['spectral_shift'] == '1':
            ctag = 'Rrs{}_{}_a'.format(qaa_wave[iw], sat)
            rrs = (coef[ctag][0]*rrs*rrs) + (coef[ctag][1]*rrs) + (coef[ctag][2])

        ## get v5/v6 indices
        if qaa_wave[iw] == 665:
            v5_idx = where(rrs < 0.0015)
            v6_idx = where(rrs >= 0.0015)

            qaa_selection = rrs*0
            if len(v5_idx[0]>0): qaa_selection[v5_idx]=5
            if len(v6_idx[0]>0): qaa_selection[v6_idx]=6
            if mask is not None: qaa_selection[where(mask)] = -1

        if mask is not None: rrs[where(mask)]=nan

        ## step 0
        ## compute Lee's subsurface rrs
        rrs = rrs/(0.52+(1.7*rrs))
        rname = 'rrs{}'.format(qaa_wave[iw])
        qaa_data[rname] = rrs
        rrs = None

        ## step 1
        ## compute Lee's u
        u = (-1.0*coef['g'][0] + sqrt((power(coef['g'][0],2.)+(4.*coef['g'][1]*qaa_data[rname]))))/(2.*coef['g'][1])
        uname = 'u{}'.format(qaa_wave[iw])
        qaa_data[uname] = u
        u = None

    ## step 2
    ## get rhow and a560 according to QAA5
    rhow_QAA5=log10((qaa_data['rrs443']+qaa_data['rrs490'])/\
                    (qaa_data['rrs560']+(5.*qaa_data['rrs665']*(qaa_data['rrs665']/qaa_data['rrs490']))))
    qaa_data['a560_QAA5'] = coef['aw'][2]+power(10.,(coef['h'][2]+(coef['h'][1]*rhow_QAA5)+coef['h'][0]*power(rhow_QAA5,2.)))

    ## get rhow and a665 according to QAA6
    rhow_QAA6=qaa_data['rrs665']/(qaa_data['rrs443']+qaa_data['rrs490'])
    qaa_data['a665_QAA6']= coef['aw'][3]+coef['k'][0]*power(rhow_QAA6,coef['k'][1])

    ## step 3
    qaa_data['bbp560_QAA5'] = ((qaa_data['u560'] * qaa_data['a560_QAA5']) / (1.0 - qaa_data['u560'])) - coef['bbw'][2]
    qaa_data['bbp665_QAA6'] = ((qaa_data['u665'] * qaa_data['a665_QAA6']) / (1.0 - qaa_data['u665'])) - coef['bbw'][3]

    ## step 4 QAAv6
    ratio = qaa_data['rrs443']/qaa_data['rrs560']
    Y = coef['l'][0] * (1.0 - coef['l'][1] * exp(coef['l'][2]*ratio))

    ## compute a QAA5 & QAA6
    for version in ['5','6']:
        ## reference wavelengths
        refw = {'5':555, '6':670}[version]
        refn = {'5':560, '6':665}[version]
        for iw,wave in enumerate(qaa_wave):
            ## get dataset tags
            bbptag = 'bbp{}_QAA{}'.format(wave, version)
            atag = 'a{}_QAA{}'.format(wave, version)
            utag = 'u{}'.format(wave)
            kdtag = 'Kd{}_QAA{}'.format(wave, version)

            ## compute bbp and q
            if not (((version == '5') & (wave == 560)) or ((version == '6') & (wave == 665))):
                ##step 5 - get bbp
                qaa_data[bbptag] = qaa_data['bbp{}_QAA{}'.format(refn,version)] * \
                                   power(float(refw)/float(wave),Y)
                ##step 6 - get a
                qaa_data[atag] = ((1.0 - qaa_data[utag]) * (coef['bbw'][iw]+qaa_data[bbptag])) / qaa_data[utag]

            ## compute Kd
            qaa_data[kdtag] = qaa_kd(qaa_data[atag], qaa_data[bbptag], coef)

    ## get switched datasets
    for iw,wave in enumerate(qaa_wave):
        for qaapar in ['a','bbp','Kd']:
            t5 = '{}{}_QAA5'.format(qaapar,wave)
            t6 = '{}{}_QAA6'.format(qaapar,wave)
            tsw = '{}{}_QAAsw'.format(qaapar,wave)
            qaa_data[tsw] = qaa_data[t6]
            if len(v5_idx[0]>0): 
                qaa_data[tsw][v5_idx]=qaa_data[t5][v5_idx]

    ## parameters for KdPAR 0-1m Nechad KdPARv2
    qaa_data['KdPAR_QAAsw']=1.2529*log(1.+qaa_data['Kd490_QAAsw'])+0.1127 

    ## get Zeu
    c1=[-0.057,0.482,4.221]
    c2=[0.183,0.702,-2.567]
    alpha=[0.090,1.465,-0.667]

    stheta = sin(ths*(pi/180))
    ctheta = cos(ths*(pi/180))
    tau = -log(0.01) # Zeu

    k1 = (c1[0] + c1[1]*sqrt(qaa_data[k_atag]) + c1[2]*qaa_data[k_bbptag])*(1+alpha[0]*stheta);
    k2 = (c2[0] + c2[1]*    (qaa_data[k_atag]) + c2[2]*qaa_data[k_bbptag])*(alpha[1]+alpha[2]*ctheta)

    # Kpar, Lee et al. 2007
    z = 1.
    qaa_data['KdPAR_QAAsw_Lee'] = k1+(k2)/sqrt(z+1.)
    qaa_data['Zeu_QAAsw_Lee'] = 4.6/qaa_data['KdPAR_QAAsw_Lee']

    ## slow (root solving)
    if compute_zeu_lee:
        y1 = (k1*k1 - k2*k2 - 2.0*tau*k1)/(k1*k1)
        y2 = (tau*tau - 2.0*tau*k1)/(k1*k1)
        y3 = (tau*tau)/(k1*k1)
        zeu = qaa_data[atag]*0.0
        zeu[:]=nan

        val = where(isfinite(y1) & isfinite(y2) & isfinite(y3))
        for i,i1 in enumerate(val[0]):
            i2=val[1][i]
            zeu[i1,i2] = real(roots((y1[i1,i2],y2[i1,i2],y3[i1,i2])))[1]
        zeu[zeu<0]=nan
        qaa_data['Zeu_QAA_Lee'] = zeu

    t1 = time.time()
    print('Finished QAA computations in {:.1f} seconds'.format(t1-t0))

    return(qaa_data)
