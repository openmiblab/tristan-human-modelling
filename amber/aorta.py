import os
import time
import pandas as pd
import numpy as np

import plot

import models.aorta_chc as ao
import dcmods.tools as tools

structure = 'aorta'

def format(xdata, ydata, weight):
    vars = {
        'weight': weight,
        'agent': 'Dotarem',
        'dose': 0.2, #mL/kg 
        'conc': 0.5, # mmol/ml 
        'rate': 1, # mL/sec
        'TR': 4.80/1000.0, #sec
        'FA': 17.0, #deg
        'R10': 1000/1463.0,
        'dt': 0.5,
    }
    pars, bounds, pfix = tools.unpack_pars(ao.PARS)
    ydata, pars, vars = ao.init(xdata, ydata, pars, vars)
    return xdata, ydata, pars, vars, pfix, bounds


def fit(time, signal, weight, path, name):

    xdata, ydata, pars, vars, pfix, bounds = format(time, signal, weight)

    # Fit model to data
    loss = tools.loss(ao.model, xdata, ydata, pars, vars=vars)
    print('Goodness of fit: ', loss)
    pars, pcov = tools.curve_fit(
        ao.model, xdata, ydata, 
        pars, vars = vars, 
        pfix = pfix,
        bounds = bounds,
        xtol = 1e-3,
    )  
    loss = tools.loss(ao.model, xdata, ydata, pars, vars=vars)
    print('Goodness of fit: ', loss)

    # Export results
    tools.to_csv(ao.PARS, pars, os.path.join(path, name+'.csv'))
    ao.plot_fit(xdata, ydata, pars, vars=vars, save=True, show=False, path=path, prefix=name)
    return ao.export_pars(pars)


def fit_data(datadir, resultsdir):
    datafile = os.path.join(datadir, 'Overview.xlsx')
    parfile = os.path.join(datadir, 'ROI size.xlsx')
    data = pd.read_excel(datafile, sheet_name='Blad1')
    const = pd.read_excel(parfile, sheet_name='Blad1')
    weight = pd.read_excel(parfile, sheet_name='Blad2')
    output = pd.DataFrame(columns=['subject','visit','structure','name','value','unit'])
    for visit in ['1','2']:
        for s in ['1','2','3','4']:
            print('Fitting ', s, visit)
            curve = 'Dog' + s + '.' + visit + ' AIF'
            wght = weight[weight.Dog==int(s)].weight.values[0]
            time = data['Time'].values
            pars = fit(time, data[curve].values, wght, resultsdir, curve)
            pars['subject'] = s
            pars['visit'] = visit
            pars['structure'] = structure
            output = pd.concat([output, pars])
    try:
        output['parameter'] = output.index
        output.to_csv(os.path.join(resultsdir, 'parameters.csv'), index=False)
        output.to_pickle(os.path.join(resultsdir, 'parameters.pkl'))
    except:
        print("Can't write to parameter file ")
        print("Please close the file before saving data")
    return os.path.join(resultsdir, 'parameters.csv')


def main(datadir, resultsdir):

    start = time.time()

    output_file = fit_data(datadir, resultsdir)

    plot.create_bar_chart(output_file, ylim={})

    print('Fit pbpk aorta 1-scan calculation time (mins): ', (time.time()-start)/60)
