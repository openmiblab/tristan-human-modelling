import os
import pandas as pd
import numpy as np

def to_tristan_units(pars):
  
    slow_time = ['BAT','BAT2','Teb','Th','Th_i','Th_f']
    for p in pars:
        if pars[p][2] == '':
            pars[p][1:] = [pars[p][1]*100, '%', pars[p][3]*100]
        if pars[p][2] == 'mL/mL':
            pars[p][1:] = [pars[p][1]*100, 'mL/100mL', pars[p][3]*100]
        if pars[p][2] == 'mL/sec/mL':
            pars[p][1:] = [pars[p][1]*6000, 'mL/min/100mL', pars[p][3]*6000]
        if p in slow_time:
            pars[p][1:] = [pars[p][1]/60, 'min', pars[p][3]/60]

    if 'CO' in pars:
        pars['CO'][1:] = [pars['CO'][1]*60/1000, 'L/min', pars['CO'][3]*60/1000]
    if 'CL' in pars:
        pars['CL'][1:] = [pars['CL'][1]*60/1000, 'L/min', pars['CL'][3]*60/1000]

    return pars

def to_df(pars):
    return pd.DataFrame.from_dict(pars, 
            orient = 'index', 
            columns = ["name", "value", "unit", 'stdev'])

def to_csv(model, file, pp=None):
    if pp is None:
        pp = model.pars()
    cols = ['symbol', "name", "initial value", 
        'value', "unit", "lower bound", "upper bound", "fix", 'stdev']
    pb = model.bounds
    mdl = model.__class__()
    perr = np.sqrt(np.diag(model.pcov))
    PARS = []
    for p, par in enumerate(model.free):
        if par in pp:
            PARS.append([par, pp[par][0], getattr(mdl,par), pp[par][1], pp[par][2], pb[0][p], pb[1][p], 0, perr[p]])
    df = pd.DataFrame(PARS, columns=cols)
    df = df.set_index('symbol')
    path = os.path.dirname(file)
    if not os.path.isdir(path):
        os.makedirs(path)
    try:
        df.to_csv(file)
    except:
        print("Can't write to file ", file)
        print("Please close the file before saving data.") 


def read_csv(file):
    try:
        df = pd.read_csv(file, index_col='symbol')   
    except:
        msg = 'Cannot read model parameters from file ' + file
        msg += '\nPlease check if the file exists and is not open in another program.'
        raise RuntimeError(msg)
    return df


def save_path(path):
    if path is None:
        path = os.path.dirname(__file__)
        path = os.path.join(path, 'results')
        if not os.path.isdir(path):
            os.mkdir(path)
    if not os.path.exists(path):
        os.makedirs(path)
    return path

