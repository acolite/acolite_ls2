## def interp2d
## interpolates 2D array for LUT lookup
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-10-30
## modifications:

def interp2d(data,xid,yid):
    dim = data.shape
    xbr = (int(xid),min(int(xid+1),dim[0]))
    ybr = (int(yid),min(int(yid+1),dim[0]))
    
    x = xid - xbr[0]
    y = yid - ybr[0]
    
    d00 = data[xbr[0],ybr[0]]
    d10 = data[xbr[1],ybr[0]]
    d01 = data[xbr[0],ybr[1]]
    d11 = data[xbr[1],ybr[1]]


    return d00 * (1-x)*(1-y) + d10 * x*(1-y) + d01 * (1-x)*y + d11 * x*y

