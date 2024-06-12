import os

import pandas as pd
import dcmri as dc

from dcmods import fig, tools

structure = 'aorta'

def fit(xdata, ydata, weight, path, name):

    # Set constants
    model = dc.AortaSignal8b()
    model.weight = weight
    model.dose = 0.2
    model.rate = 1
    model.R10 = 1000/1463.0
    model.TR = 4.80/1000.0
    model.FA = 17.0
    model.dose_tolerance = 0.1
    model.agent = 'Dotarem'
    model.dt = 0.5

    # Get initial values
    model.initialize('TRISTAN')
    model.pretrain(xdata, ydata)
    loss = model.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss)

    # Fit model to data
    model.train(xdata, ydata, bounds='TRISTAN', xtol=1e-3)
    loss = model.cost(xdata, ydata)
    print('Goodness of fit (optimal): ', loss)

    # Plot fits
    kwargs = {'color':['lightcoral','darkred'], 'label':'Aorta', 
              'show':False, 'save':True, 'path':path, 'prefix':name}
    fig.tissue(model, xdata, ydata, **kwargs)
    BAT = model.pars[0]
    fig.tissue(model, xdata, ydata, win='pass1', xlim=[BAT-20, BAT+160], **kwargs)
    
    # Export results
    tools.to_csv(model, os.path.join(path, name + '.csv'))
    pars = model.pfree(units='custom') + model.pdep(units='custom')
    return tools.to_df(pars)


def main(datadir, resultsdir):

    datafile = os.path.join(datadir, 'Overview.xlsx')
    parfile = os.path.join(datadir, 'ROI size.xlsx')
    data = pd.read_excel(datafile, sheet_name='Blad1')
    const = pd.read_excel(parfile, sheet_name='Blad1')
    weight = pd.read_excel(parfile, sheet_name='Blad2')
    output = None
    for visit in ['1','2']:
        for s in ['1','2','3','4','5','6','7','8']:
            print('Fitting ', s, visit)
            curve = 'Dog' + s + '.' + visit + ' AIF'
            wght = weight[weight.Dog==int(s)].weight.values[0]
            time = data['Time'].values
            pars = fit(time, data[curve].values, wght, resultsdir, curve)
            pars['subject'] = s
            pars['visit'] = visit
            pars['structure'] = structure
            if output is None:
                output = pars
            else:
                output = pd.concat([output, pars])
    try:
        output = output.reindex(columns=['subject','visit','structure','name','value','unit'])
        output['parameter'] = output.index
        output.to_csv(os.path.join(resultsdir, 'parameters.csv'), index=False)
        output.to_pickle(os.path.join(resultsdir, 'parameters.pkl'))
    except:
        print("Can't write to parameter file ")
        print("Please close the file before saving data")

