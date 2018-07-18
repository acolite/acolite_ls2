## def metadata_read
## reads landsat metadata
## written by Quinten Vanhellemont, RBINS
## 2017-04-13
## modifications:

def metadata_read(metafile):
    mdata={}
    with open(metafile,'r') as f:
        for line in f.readlines():
            split = line.split('=')
            for s, sp in enumerate(split):
                split[s] = sp.strip(' \n "')
            if len(split) is not 2: continue
                
            if split[0] == 'GROUP':
                cur_group = split[1]
                #if cur_group == 'L1_METADATA_FILE': continue
                group_data = {}
            elif split[0] == 'END_GROUP':
                #if split[1] == 'L1_METADATA_FILE': continue
                mdata[cur_group] = group_data
            else:
                group_data[split[0]] = split[1]
    return(mdata)
