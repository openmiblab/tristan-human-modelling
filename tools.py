import os
import time
import pickle

import pandas as pd
import numpy as np

import plot

AORTA_PARS = ['BAT','CO','Thl','Dhl','To','Eb','Eo','Toe','BAT2', 
              'AUC_R1b','AUC_Cb', 'AUC35_R1b','AUC35_Cb']

def compute(fit_subj, resultspath):

    start = time.time()
    resultspath = save_path(resultspath)

    with open(os.path.join(resultspath, 'data.pkl'), 'rb') as fp:
        data_dict = pickle.load(fp)

    output = None
    for visit in data_dict:
        for subj in data_dict[visit]:
            print('Fitting aorta and liver of ', visit, subj)
            name = subj + '_' + visit
            pars = fit_subj(data_dict[visit][subj], resultspath, name)
            pars['subject'] = subj
            pars['visit'] = visit
            structure = []
            for p in pars.index.values:
                if p in AORTA_PARS:
                    structure.append('aorta')
                else:
                    structure.append('liver')
            pars['structure'] = structure
            if output is None:
                output = pars
            else:
                output = pd.concat([output, pars])

    # Format output and save
    cols = ['subject','visit','structure','name','value','unit','stdev']
    output = output.reindex(columns=cols)
    output['parameter'] = output.index
    output.to_csv(os.path.join(resultspath, 'parameters.csv'), index=False)
    output.to_pickle(os.path.join(resultspath, 'parameters.pkl'))
    
    print('Calculation time (mins): ', (time.time()-start)/60)


def compute_vart(fit_subj, resultspath, datapath, 
                 acq_times = [5,10,15,20,25,30,35,40]):

    start = time.time()
    resultspath = save_path(resultspath)

    with open(os.path.join(datapath, 'data.pkl'), 'rb') as fp:
        data_dict = pickle.load(fp)

    output = None
    for visit in data_dict:
        for subj in data_dict[visit]:
            for tacq in acq_times:
                print('Fitting data of ', visit, subj, tacq)
                name = subj + '_' + visit + '_' + str(tacq)
                pars = fit_subj(
                    data_dict[visit][subj], resultspath, name, tacq*60)
                pars['subject'] = subj
                pars['visit'] = visit
                structure = []
                for p in pars.index.values:
                    if p in AORTA_PARS:
                        structure.append('aorta')
                    else:
                        structure.append('liver')
                pars['structure'] = structure
                pars['tacq'] = tacq
                if output is None:
                    output = pars
                else:
                    output = pd.concat([output, pars])

    # Format output and save
    output = output.reindex(columns= ['subject','visit','structure','name',
                                      'value','unit','stdev','tacq'])
    output['parameter'] = output.index
    output.to_csv(os.path.join(resultspath, 'parameters.csv'), index=False)
    output.to_pickle(os.path.join(resultspath, 'parameters.pkl'))
    
    print('Calculation time (mins): ', (time.time()-start)/60)


def to_tristan_units(pars):
  
    slow_time = ['BAT','BAT2','Toe','Th','Th_i','Th_f']
    for p in pars:
        if pars[p][2] == '':
            pars[p][1:] = [pars[p][1]*100, '%', pars[p][3]*100]
        if pars[p][2] == 'mL/cm3':
            pars[p][1:] = [pars[p][1]*100, 'mL/100cm3', pars[p][3]*100]
        if pars[p][2] == 'mL/sec/cm3':
            pars[p][1:] = [pars[p][1]*6000, 'mL/min/100cm3', pars[p][3]*6000]
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
        pp = model.export_params()
    cols = ['symbol', "name", "initial value", 
        'value', "unit", "lower bound", "upper bound", "fix", 'stdev']
    perr = np.sqrt(np.diag(model.pcov))
    PARS = []
    for p, par in enumerate(model.free):
        if par in pp:
            PARS.append([
                par, 
                pp[par][0], 
                model._params()[par]['init'],
                pp[par][1], 
                pp[par][2], 
                model.free[par][0], 
                model.free[par][1], 
                0, 
                perr[p]])
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
        raise RuntimeError(
            'Cannot read model parameters from file ' + file + '\nPlease check'
            ' if the file exists and is not open in another program.')
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

