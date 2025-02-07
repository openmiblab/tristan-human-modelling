import os
import pickle

import numpy as np
import dcmri as dc

import tools
from tristan import data


def params(model:dc.AortaLiver, tb, Sb, tl, Sl):

     # Compute AUC over 3hrs
    model.tmax = model.BAT+180*60
    t, cb, Cl = model.conc()
    t, R1b, R1l = model.relax()
    AUC_DR1b = np.trapezoid(R1b-model.R10a, t)
    AUC_Cb = np.trapezoid(cb, model.t) 
    AUC_DR1l = np.trapezoid(R1l-model.R10l, t)
    AUC_Cl = np.trapezoid(Cl, t)

    # Compute relative enhancement at 20mins
    tRE = model.BAT + 20*60
    RE_R1b = (R1b[t<tRE][-1] - R1b[0])/R1b[0]
    RE_R1l = (R1l[t<tRE][-1] - R1l[0])/R1l[0]
    S0b = np.mean(Sb[tb<model.BAT-30])
    S0l = np.mean(Sl[tl<model.BAT-30])
    RE_Sb = (Sb[tb<tRE][-1] - S0b)/S0b
    RE_Sl = (Sl[tl<tRE][-1] - S0l)/S0l


     # Compute AUC over 35min
    model.tmax = model.BAT+35*60
    t, cb, Cl = model.conc()
    t, R1b, R1l = model.relax()
    AUC35_DR1b = np.trapezoid(R1b-model.R10a, t)
    AUC35_Cb = np.trapezoid(cb, model.t) 
    AUC35_DR1l = np.trapezoid(R1l-model.R10l, t)
    AUC35_Cl = np.trapezoid(Cl, t) 

    pars = model.export_params()
    pars['AUC_R1b']=['AUC for DR1b (0-inf)', AUC_DR1b, '',0]
    pars['AUC_Cb']=['AUC for Cb (0-inf)', 1000*AUC_Cb, 'mM*sec',0]
    pars['AUC_R1l']=['AUC for DR1l (0-inf)', AUC_DR1l, '',0]
    pars['AUC_Cl']=['AUC for Cl (0-inf)', 1000*AUC_Cl, 'mM*sec',0]  
    pars['AUC35_R1b']=['AUC for DR1b (0-35min)', AUC35_DR1b, '',0]
    pars['AUC35_Cb']=['AUC for Cb (0-35min)', 1000*AUC35_Cb, 'mM*sec',0]
    pars['AUC35_R1l']=['AUC for DR1l (0-35min)', AUC35_DR1l, '',0]
    pars['AUC35_Cl']=['AUC for Cl (0-35min)', 1000*AUC35_Cl, 'mM*sec',0]  
    pars['RE_R1b']=['RE for R1b at 20min', 100*RE_R1b, '%',0]
    pars['RE_R1l']=['RE for R1l at 20min', 100*RE_R1l, '%',0]
    pars['RE_Sb']=['RE for Sb at 20min', 100*RE_Sb, '%',0]
    pars['RE_Sl']=['RE for Sl at 20min', 100*RE_Sl, '%',0]       

    return pars  


def figure(model:dc.AortaLiver, 
            xdata:tuple[np.ndarray, np.ndarray], 
            ydata:tuple[np.ndarray, np.ndarray], 
            path, name, t, R1a, R1l):
    file = os.path.join(path, name)
    ya = [dc.signal_ss(model.S0a, R1a[0], model.TR, model.FA),
          dc.signal_ss(model.S0a, R1a[1], model.TR, model.FA)]
    yl = [dc.signal_ss(model.S0l, R1l[0], model.TR, model.FA),
          dc.signal_ss(model.S0l, R1l[1], model.TR, model.FA)]
    test=((t,ya),(t,yl))
    BAT = model.BAT
    model.plot(xdata, ydata, 
               fname=file + '.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+1200], 
               fname=file + '_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_win3.png', ref=test, show=False) 



def fit_subj(data, path, name, tacq=None):

    xdata = data['xdata']
    ydata = data['ydata']

    # Truncate data if requested
    if tacq is not None:
        idx0, idx1 = xdata[0]<tacq, xdata[1]<tacq
        xdata = (xdata[0][idx0], xdata[1][idx1])
        ydata = (ydata[0][idx0], ydata[1][idx1])
    
    # Fit model to data
    model = dc.AortaLiver(**data['params'])
    loss0 = model.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss0)
    model.train(xdata, ydata, xtol=1e-3, verbose=2)
    loss1 = model.cost(xdata, ydata)
    print('Goodness of fit (improvement, %): ', 100*(loss0-loss1)/loss0)

    # Export data
    figure(model, xdata, ydata, path, name, 
           data['tR1'], data['R1a'], data['R1l'])
    pars = params(model, xdata[0], ydata[0], xdata[1], ydata[1])
    tools.to_csv(model, os.path.join(path, name + '.csv'), pars)
    pars = tools.to_tristan_units(pars)
    return tools.to_df(pars)



def format_data(datapath, resultspath):

    resultspath = tools.save_path(resultspath)

    data_dict = {}
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        data_dict[visit] = {}
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            subj_data = data.read(subj)
            data_dict[visit][s[:3]] = {
                'xdata': (subj_data['xdata'][0], subj_data['xdata'][2]),
                'ydata': (subj_data['ydata'][0], subj_data['ydata'][2]),
                'tR1':  subj_data['tR1'][:2],
                'R1a':  subj_data['R1a'][:2],
                'R1l':  subj_data['R1l'][:2],
                'params': {
                    'weight':  subj_data['weight'],
                    'agent': 'gadoxetate',
                    'dose':  subj_data['dose1'],
                    'rate':  1,
                    'field_strength':  3.0,
                    't0':  subj_data['baseline'],
                    'TR':  3.71/1000.0,
                    'FA':  15,
                    'TS':  subj_data['time1'][1]-subj_data['time1'][0],
                    'R10a': subj_data['R1a'][0],
                    'R10l': subj_data['R1l'][0],
                    'H':  0.45,
                    'vol':  subj_data['liver_volume'],
                }
            }

    with open(os.path.join(resultspath, 'data.pkl'), 'wb') as fp:
        pickle.dump(data_dict, fp, protocol=pickle.HIGHEST_PROTOCOL)

