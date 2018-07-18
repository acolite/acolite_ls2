## def interp3d
## interpolates 3D array for LUT lookup
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-05
## modifications:

def interp3d(data,xid,yid,zid):
    dim = data.shape
    xbr = (int(xid),min(int(xid+1),dim[0]))
    ybr = (int(yid),min(int(yid+1),dim[0]))
    zbr = (int(zid),min(int(zid+1),dim[0]))
    
    x = xid - xbr[0]
    y = yid - ybr[0]
    z = zid - zbr[0]
    
    d000 = data[xbr[0],ybr[0],zbr[0]]
    d100 = data[xbr[1],ybr[0],zbr[0]]
    d010 = data[xbr[0],ybr[1],zbr[0]]
    d001 = data[xbr[0],ybr[0],zbr[1]]
    d101 = data[xbr[1],ybr[0],zbr[1]]
    d011 = data[xbr[0],ybr[1],zbr[1]]
    d110 = data[xbr[1],ybr[1],zbr[0]]
    d111 = data[xbr[1],ybr[1],zbr[1]]

    return d000 * (1-x)*(1-y)*(1-z) + d100 * x*(1-y)*(1-z) + d010 * (1-x)*y*(1-z) + d001 * (1-x)*(1-y)*z + d101 * x*(1-y)*z + d011 * (1-x)*y*z + d110 * x*y*(1-z) + d111 * x*y*z

