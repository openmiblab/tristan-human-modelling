import os
import time

import numpy as np
import pandas as pd

import dcmri

import tristan.data as data
import plot

import dcmods.tools as tools
import models.liver_2cm_2scan as liv_2cm
import models.liver_2cm_ns_2scan as liv_2cm_ns
import models.aorta_pcc_2scan as ao
import tristan.aorta_2scan as aofit

def format(data, aortapars):

    # Get aorta data
    df = tools.read_csv(aortapars)
    pars = df.value.values
    xdata, _, _, vars, _, _, _ = aofit.format(data)
    t, cb = ao.model(xdata, *pars, return_conc=True, **vars)
    BAT1, BAT2 = pars[0], pars[1]

    # Format liver data
    liv = liv_2cm if data['dose2']==0 else liv_2cm_ns
    xvalid = np.concatenate([data['liver_valid1'], data['liver_valid2']])
    xdata = np.concatenate([data['time1'], data['time2']])
    ydata = np.concatenate([data['liver1'], data['liver2']])
    pars, bounds, pfix = tools.unpack_pars(liv.PARS)
    vars = {
        'dt': t[1]-t[0],
        'cb': cb,
        'R10': 1000.0/data['T1liver1'],
        'TR': 3.71/1000.0,
        'FA': 15.0,
        'tR12': data['T1time3'],
    }
    pars, vars = liv.init(xdata, ydata, pars, vars, BAT1, 1000.0/data['T1liver3'])
    return xdata, ydata, pars, vars, pfix, xvalid, bounds, BAT1, BAT2, liv


def fit(data, path, name, aortapars):

    xdata, ydata, pars, vars, pfix, xvalid, bounds, BAT1, BAT2, liv = format(data, aortapars)

    # Fit model to data
    loss = tools.loss(liv.model, xdata, ydata, 
        pars, vars=vars, xvalid=xvalid)
    print('Goodness of fit: ', loss)
    pars, pcov = tools.curve_fit(
        liv.model, xdata, ydata, 
        pars, vars = vars, 
        pfix = pfix,
        xvalid = xvalid,
        bounds = bounds,
        xtol = 1e-6,
    )  
    loss = tools.loss(liv.model, xdata, ydata,
        pars, vars=vars, xvalid=xvalid)
    print('Goodness of fit: ', loss)

    # Export results
    xcheck = np.array([0, data['T1time2'], data['T1time3']])
    R1 = [1000.0/data['T1liver1'], 1000.0/data['T1liver2'], 1000.0/data['T1liver3']]
    ycheck = [
        dcmri.signal_spgress(vars['TR'], vars['FA'], R1[0], pars[-2]),
        dcmri.signal_spgress(vars['TR'], vars['FA'], R1[1], pars[-2]),
        dcmri.signal_spgress(vars['TR'], vars['FA'], R1[2], pars[-1]),
    ]
    tools.to_csv(liv.PARS, pars, os.path.join(path, name + '.csv'))
    liv.plot_fit(xdata, ydata, pars, BAT1, BAT2, vars=vars, xvalid=xvalid, save=True, show=False, path=path, prefix=name+'_Liver', xcheck=xcheck, ycheck=ycheck)
    if data['dose2']==0:
        return liv.export_pars(pars, data['liver_volume'])
    else:
        return liv.export_pars(pars, data['liver_volume'], vars['dt'], data['t0'], xdata, vars['tR12'])


def main(datapath, results):

    start = time.time()
    resultspath = os.path.join(results, 'liver_2scan')
    aortaresults = os.path.join(results, 'aorta_2scan') 
    output_file = os.path.join(resultspath, 'parameters.csv')

    output = None
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            print('Fitting liver of ', visit, subj)
            subj_data = data.read(subj)
            name = s[:3] + '_' + visit
            aortapars = os.path.join(aortaresults, name + '.csv')
            liver_pars = fit(subj_data, resultspath, name, aortapars)
            liver_pars['subject'] = s[:3]
            liver_pars['visit'] = visit
            liver_pars['structure'] = 'liver'
            if output is None:
                output = liver_pars
            else:
                output = pd.concat([output, liver_pars])
    try:
        output = output.reindex(columns=['subject','visit','structure','name','value','unit'])
        output['parameter'] = output.index
        output.to_csv(output_file, index=False)
        plot.create_bar_chart(output_file)
    except:
        print("Can't write to file ", output_file)
        print("Please close the file before saving data")
    
    print('Calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()