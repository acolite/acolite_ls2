## def widx
##
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-01-30
## modifications:

def widx(waves, wv):
    return(min(enumerate(waves), key=lambda x: abs(x[1]-wv)))
