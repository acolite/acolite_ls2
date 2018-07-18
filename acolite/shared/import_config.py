## simple function to import txt config file
## QV 2017-05-24

def import_config(file):
    config={}
    with open(file, 'r') as f:
        for line in f.readlines():
            if line[0] == '#': continue
            split=line.strip('\n').split('=')
            if len(split) != 2: continue
            config[split[0]]=split[1]
    return(config)
