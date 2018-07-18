## def similarity_ratio_wave
## returns similarity spectrum ratio for given wavelengths
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-04-17
## modifications: 
##                2018-07-18 (QV) changed acolite import name

def similarity_ratio_wave(w1, w2):
    import acolite as ac
    
    ## get similarity spectrum table
    ssd = ac.shared.similarity_read()
    idr, wr = ac.shared.closest_idx(ssd['wave'], w1/1000.)
    idn, wn = ac.shared.closest_idx(ssd['wave'], w2/1000.)
    return(ssd['ave'][idr]/ssd['ave'][idn])
