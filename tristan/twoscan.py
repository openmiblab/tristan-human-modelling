import os
import time

import numpy as np
import pandas as pd
import dcmri as dc

import plot, tools
from tristan import data

def params(model:dc.AortaLiver2scan, xdata, t0):
    model.dose[1] = 0

    model.tmax = model.BAT+180*60
    t, cb, Cl = model.conc()
    t, R1b, R1l = model.relax()
    AUC_DR1b = np.trapz(R1b-model.R10b, t)
    AUC_Cb = np.trapz(cb, model.t) 
    AUC_DR1l = np.trapz(R1l-model.R10l, t)
    AUC_Cl = np.trapz(Cl, t)

    model.tmax = model.BAT+35*60
    t, cb, Cl = model.conc()
    t, R1b, R1l = model.relax()
    AUC35_DR1b = np.trapz(R1b-model.R10b, t)
    AUC35_Cb = np.trapz(cb, model.t) 
    AUC35_DR1l = np.trapz(R1l-model.R10l, t)
    AUC35_Cl = np.trapz(Cl, t) 

    pars = model.pars()
    pars['AUC_R1b']=['AUC for DR1b (0-inf)', AUC_DR1b, '',0]
    pars['AUC_Cb']=['AUC for Cb (0-inf)', 1000*AUC_Cb, 'mM*sec',0]
    pars['AUC_R1l']=['AUC for DR1l (0-inf)', AUC_DR1l, '',0]
    pars['AUC_Cl']=['AUC for Cl (0-inf)', 1000*AUC_Cl, 'mM*sec',0]  
    pars['AUC35_R1b']=['AUC for DR1b (0-35min)', AUC35_DR1b, '',0]
    pars['AUC35_Cb']=['AUC for Cb (0-35min)', 1000*AUC35_Cb, 'mM*sec',0]
    pars['AUC35_R1l']=['AUC for DR1l (0-35min)', AUC35_DR1l, '',0]
    pars['AUC35_Cl']=['AUC for Cl (0-35min)', 1000*AUC35_Cl, 'mM*sec',0]       

    # Timings needed for plotting etc
    pars['t0']=["Start time first acquisition", (t0+xdata[0][0])/(60*60), 'hrs',0]
    pars['t1']=["End time first acquisition", (t0+xdata[0][-1])/(60*60), 'hrs',0]
    pars['t2']=["Start time second acquisition", (t0+xdata[1][0])/(60*60), 'hrs',0]
    pars['t3']=["End time second acquisition", (t0+xdata[1][-1])/(60*60), 'hrs',0]
    pars['dt1']=["Time step first acquisition", model.tacq, 'sec',0]
    pars['dt2']=["Time step second acquisition", model.tacq2, 'sec',0]

    return pars  


def figure(model:dc.AortaLiver2scan, 
            xdata:tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], 
            ydata:tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray], 
            path, name, t, R1a, R1l):
    ya = [dc.signal_ss(R1a[0], model.S0b, model.TR, model.FA),
          dc.signal_ss(R1a[1], model.S0b, model.TR, model.FA),
          dc.signal_ss(R1a[2], model.S02b, model.TR, model.FA)]
    yl = [dc.signal_ss(R1l[0], model.S0l, model.TR, model.FA),
          dc.signal_ss(R1l[1], model.S0l, model.TR, model.FA),
          dc.signal_ss(R1l[2], model.S02l, model.TR, model.FA)]
    test = ((t,ya),(t,yl))
    file = os.path.join(path, name)
    model.plot(xdata, ydata, fname=file + '__win0.png', testdata=test, show=False)
    BAT = model.BAT
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+1200], fname=file + '_scan1_win1.png', testdata=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], fname=file + '_scan1_win2.png', testdata=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], fname=file + '_scan1_win3.png', testdata=test, show=False)
    BAT = model.BAT2
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+1200], fname=file + '_scan2_win1.png', testdata=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], fname=file + '_scan2_win2.png', testdata=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], fname=file + '_scan2_win3.png', testdata=test, show=False)


def fit(data, path, name):

    xdata, ydata = data['xdata'], data['ydata']

    # Get model
    model = dc.AortaLiver2scan(
        tacq = data['time1'][1]-data['time1'][0],
        tacq2 = data['time2'][1]-data['time2'][0],
        weight = data['weight'],
        agent = 'gadoxetate',
        dose = [data['dose1'], data['dose2']],
        rate = 1,
        field_strength = 3.0,
        t0 = data['baseline'],
        TR = 3.71/1000.0,
        FA = 15,
        R10b = data['R1a'][0],
        R10l = data['R1l'][0],
        R102b = data['R1a'][2],
        R102l = data['R1l'][2],
        Hct = 0.45,
        vol = data['liver_volume'],
    )

    # Fit model to data
    loss0 = model.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss0)
    model.train(xdata, ydata, xtol=1e-3, verbose=2)
    loss1 = model.cost(xdata, ydata)
    print('Goodness of fit (improvement, %): ', 100*(loss0-loss1)/loss0)

    # Export data
    figure(model, xdata, ydata, path, name, 
           data['tR1'], data['R1a'], data['R1l'])
    pars = params(model, (xdata[0],xdata[1]), data['t0'])
    tools.to_csv(model, os.path.join(path, name + '.csv'), pars)
    pars = tools.to_tristan_units(pars)
    return tools.to_df(pars)


def main(datapath, results):

    start = time.time()
    resultspath = os.path.join(results, 'twoscan') 
    resultspath = tools.save_path(resultspath)

    output = None
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            print('Fitting aorta and liver of ', visit, subj)
            subj_data = data.read(subj)
            name = s[:3] + '_' + visit
            pars = fit(subj_data, resultspath, name)
            pars['subject'] = s[:3]
            pars['visit'] = visit
            structure = []
            for p in pars.index.values:
                if p in ['BAT','CO','Thl','Dhl','To',
                         'Eb','Eo','Teb','BAT2','Tc',
                         'AUC_R1b','AUC_Cb','AUC35_R1b','AUC35_Cb']:
                    structure.append('aorta')
                else:
                    structure.append('liver')
            pars['structure'] = structure
            if output is None:
                output = pars
            else:
                output = pd.concat([output, pars])

    # Format output and save
    output = output.reindex(columns=['subject','visit','structure','name','value','unit','stdev'])
    output['parameter'] = output.index
    output.to_csv(os.path.join(resultspath, 'parameters.csv'), index=False)
    output.to_pickle(os.path.join(resultspath, 'parameters.pkl'))
    plot.create_bar_chart(os.path.join(resultspath, 'parameters.pkl'))
    
    print('Calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()