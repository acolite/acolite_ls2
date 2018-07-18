## def rmsd
## calculates rmsd
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2016-07-11
## modifications:

def rmsd(x, y):
    if len(x) != len(y): return 1
    n = len(x)
    data = [pow((x[i]-y[i]),2) for i,n in enumerate(x)]
    return pow(sum(data)/float(n),0.5)

