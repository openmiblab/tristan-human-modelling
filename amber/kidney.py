import os

import pandas as pd
import numpy as np
import dcmri as dc

from dcmods import tools, fig


def fit(xdata, ydata, weight, R1, kidney_volume, aortapars, path, name):

    # Get aorta data
    df = tools.read_csv(aortapars)
    dt = 0.5
    aorta = dc.AortaSignal8b()
    aorta.pars = df.value.values
    aorta.weight = weight
    aorta.dose = 0.2
    aorta.rate = 1
    aorta.dose_tolerance = 0.1
    aorta.agent = 'Dotarem'
    aorta.dt = dt
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    cb = aorta.predict(t, return_conc=True)
    BAT = aorta.pars[0]

    # For kidney model
    kidney = dc.KidneySignal9()
    kidney.dt = dt
    kidney.J_aorta = cb*aorta.pars[1]
    kidney.TR = 4.80/1000.0
    kidney.FA = 17
    kidney.R10 = R1
    kidney.kidney_volume = kidney_volume
    kidney.BAT = BAT
    kidney.CO = aorta.pars[1]
    kidney.Hct = 0.36

    # Fit model to data
    kidney.initialize('Amber')
    kidney.pretrain(xdata, ydata)
    loss = kidney.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss)

    kidney.train(xdata, ydata, bounds='Amber', xtol=1e-3)
    loss = kidney.cost(xdata, ydata)
    print('Goodness of fit (optimal): ', loss)

    # Export results
    kwargs = {'color':['cornflowerblue','darkblue'], 'label':['Plasma', 'Tubuli'], 
              'show':False, 'save':True, 'path':path, 'prefix':name+'_kidney'}
    fig.tissue_2cm(kidney, xdata, ydata, **kwargs)
    fig.tissue_2cm(kidney, xdata, ydata, win='win', xlim=[BAT-20, BAT+40], **kwargs)
    fig.tissue_2cm(kidney, xdata, ydata, win='win_', xlim=[BAT-20, BAT+160], **kwargs)
    
    # Export data
    tools.to_csv(kidney, os.path.join(path, name + '.csv'), 'Amber')
    pars = kidney.pfree(units='custom') + kidney.pdep(units='custom')
    return tools.to_df(pars)


def main(datadir, aorta_results, kidney_results):
    
    datafile = os.path.join(datadir, 'Overview.xlsx')
    parfile = os.path.join(datadir, 'ROI size.xlsx')
    data = pd.read_excel(datafile, sheet_name='Blad1')
    const = pd.read_excel(parfile, sheet_name='Blad1')
    weight = pd.read_excel(parfile, sheet_name='Blad2')
    output = None
    for s in ['1','2','3','4','5','6','7','8']:
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

