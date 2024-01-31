import os
import time
import pandas as pd
import numpy as np

import dcmri

import tristan.data as data
import plot

import models.aorta_pcc as ao
import dcmods.tools as tools


def format(data):
    pars, bounds, pfix = tools.unpack_pars(ao.PARS)
    vars = {
        'weight': data['weight'],
        'dose': data['dose1'],
        'rate': 1, # mL/sec
        'R10': 1000.0/data['T1aorta1'],
        'TR': 3.71/1000.0,
        'FA': 15.0,
        'dose_tolerance': 0.1,
        'dt': 0.5,
    }
    xdata = np.array(data['time1'])
    ydata = np.array(data['aorta1'])
    xvalid = np.array(data['aorta_valid1'])
    pars, vars = ao.init(xdata, ydata, pars, vars)
    return xdata, ydata, pars, vars, pfix, xvalid, bounds


def fit(data, path, name):

    # Format data
    xdata, ydata, pars, vars, pfix, xvalid, bounds = format(data)

    # Fit model to data
    loss = tools.loss(ao.model, xdata, ydata, pars, vars=vars, xvalid=xvalid)
    print('Goodness of fit: ', loss)
    pars, pcov = tools.curve_fit(
        ao.model, xdata, ydata, 
        pars, vars = vars, 
        pfix = pfix,
        xvalid = xvalid,
        bounds = bounds,
        xtol = 1e-3,
    )  
    loss = tools.loss(ao.model, xdata, ydata, pars, vars=vars, xvalid=xvalid)
    print('Goodness of fit: ', loss)

    # Export results
    xcheck = np.array([0, data['T1time2']])
    R1 = [1000.0/data['T1aorta1'], 1000.0/data['T1aorta2']]
    ycheck = [
        dcmri.signal_spgress(vars['TR'], vars['FA'], R1[0], vars['S0']),
        dcmri.signal_spgress(vars['TR'], vars['FA'], R1[1], vars['S0']),
    ]
    tools.to_csv(ao.PARS, pars, os.path.join(path, name + '.csv'))
    ao.plot_fit(xdata, ydata, pars, vars=vars, xvalid=xvalid, xcheck=xcheck, ycheck=ycheck, save=True, show=False, path=path, prefix=name)
    return ao.export_pars(pars)


def main(datapath, results):

    start = time.time()
    resultspath = os.path.join(results, 'aorta_1scan')
    output_file = os.path.join(resultspath, 'parameters.csv')

    output = None
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            print('Fitting aorta of ', visit, subj)
            subj_data = data.read(subj)
            name = s[:3] + '_' + visit
            aorta_pars = fit(subj_data, resultspath, name)
            aorta_pars['subject'] = s[:3]
            aorta_pars['visit'] = visit
            aorta_pars['structure'] = 'aorta'
            if output is None:
                output = aorta_pars
            else:
                output = pd.concat([output, aorta_pars])
    try:
        output = output.reindex(columns=['subject','visit','structure','name','value','unit'])
        output['parameter'] = output.index
        output.to_csv(output_file, index=False)
        plot.create_bar_chart(output_file, ylim={})
    except:
        print("Can't write to file ", output_file)
        print("Please close the file before saving data")

    print('Fit aorta 1-scan calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()