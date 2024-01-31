import os
import time
import pandas as pd
import numpy as np

import plot

import dcmods.tools as tools
import models.kidney_nephron_short_pbpk as kid
import models.aorta_chc as ao
import amber.aorta as aofit


def format(time, signal, weight, R1, kidney_volume, aortapars):

    # Get aorta data
    df = tools.read_csv(aortapars)
    pars = df.value.values
    _, _, _, aovars, _, _ = aofit.format(time, signal, weight)
    t, cb = ao.model(time, *pars, return_conc=True, **aovars)
    BAT = pars[0]
    CO = pars[1]

    # Format kidney data
    vars = {
        'R10k': R1,
        'dt': aovars['dt'],
        'J_aorta': CO*cb/1000,
        'TR': aovars['TR'],
        'FA': aovars['FA'],
        'agent': aovars['agent'],
        'kidney_volume': kidney_volume,
    }
    pars, bounds, pfix = tools.unpack_pars(kid.PARS)
    xdata, ydata = time, signal
    vars = kid.init(xdata, ydata, vars, BAT)
    return xdata, ydata, pars, vars, pfix, bounds, BAT, CO


def fit(time, signal, weight, R1, kidney_volume, aortapars, path, name):

    xdata, ydata, pars, vars, pfix, bounds, BAT, CO = format(time, signal, weight, R1, kidney_volume, aortapars)
    
    # Fit model to data
    loss = tools.loss(kid.model, xdata, ydata, pars, vars=vars)
    print('Goodness of fit: ', loss)
    pars, pcov = tools.curve_fit(
        kid.model, xdata, ydata, 
        pars, vars = vars, 
        pfix = pfix,
        bounds = bounds,
        xtol = 1e-6,
    )  
    loss = tools.loss(kid.model, xdata, ydata, pars, vars=vars)
    print('Goodness of fit: ', loss)

    # Export results
    Hct = 0.36
    tools.to_csv(kid.PARS, pars, os.path.join(path, name + '.csv'))
    kid.plot_fit(xdata, ydata, pars, BAT, vars=vars, save=True, show=False, path=path, prefix=name)
    return kid.export_pars(pars, kidney_volume, time.max(), CO, Hct)


def fit_data(datadir, aorta_results, kidney_results):
    datafile = os.path.join(datadir, 'Overview.xlsx')
    parfile = os.path.join(datadir, 'ROI size.xlsx')
    data = pd.read_excel(datafile, sheet_name='Blad1')
    const = pd.read_excel(parfile, sheet_name='Blad1')
    weight = pd.read_excel(parfile, sheet_name='Blad2')
    output = None
    for s in ['1','2','3','4']:
        for visit in ['1','2']:
            for kid in ['LK','RK']:
                print('Fitting ', s, visit, kid)
                curve = 'Dog' + s + '.' + visit + ' ' + kid
                wght = weight[weight.Dog==int(s)].weight.values[0]
                time = data['Time'].values
                aortapars = os.path.join(aorta_results, 'Dog' + s + '.' + visit + ' AIF.csv')
                c = const[(const.Dog==int(s)) & (const.Session==int(visit)) & (const.Kidney==kid)]
                kidney_volume = c['ROI size (ml)'].values[0]
                R1 = 1000/c['T1 Kidney'].values[0]
                pars = fit(time, data[curve].values, wght, R1, kidney_volume, aortapars, kidney_results, curve)
                pars['subject'] = s[:3]
                pars['visit'] = visit
                pars['structure'] = kid
                if output is None:
                    output = pars
                else:
                    output = pd.concat([output, pars])
    try:
        output = output.reindex(columns=['subject','visit','structure','name','value','unit'])
        output['parameter'] = output.index
        output_file = os.path.join(kidney_results, 'parameters.csv')
        output.to_csv(output_file, index=False)
        output.to_pickle(os.path.join(kidney_results, 'parameters.pkl'))
    except:
        print("Can't write to parameter file")
        print("Please close the file before saving data")

    return output_file



def main(datadir, aorta_results, kidney_results):
    start = time.time()
    output_file = fit_data(datadir, aorta_results, kidney_results)
    plot.create_bar_chart(output_file, ylim={})
    print('Calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()