## def lutpos
## finds position of value in a (sorted) vector for LUT lookup
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:

def lutpos(vector, value):
    from numpy import searchsorted
    uidx = searchsorted(vector, value, side='right')
    uidx = max(0,uidx)
    uidx = min(uidx,len(vector)-1)
    uvalue = vector[uidx]
    
    if uvalue > value: 
        lidx = uidx-1
        index = lidx + (value - (vector[lidx])) / abs(vector[lidx]-vector[uidx])
        bracket=(lidx,uidx)
        return index, bracket
    if uvalue <= value:
        lidx = uidx
        uidx = uidx
        index = float(lidx)
        bracket=(lidx,uidx)
        return index, bracket
    
    return 0, (0,0)
