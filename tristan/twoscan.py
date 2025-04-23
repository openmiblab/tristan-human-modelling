import os
import time

import numpy as np
import dcmri as dc

from tristan import tools



def compute(datapath, resultspath):

    start = time.time()

    if not os.path.exists(resultspath):
        os.makedirs(resultspath)

    rois, pars = dc.read_dmr(datapath, nest=True, valsonly=True)
    for subj in rois.keys():
        for visit in rois[subj].keys():
            name = subj + '_' + visit
            fit_subj(rois[subj][visit], pars[subj][visit], resultspath, name)

    print('Calculation time (mins): ', (time.time()-start)/60)


def fit_subj(rois, pars, path, name):

    aorta_1 = rois['aorta_1_accept']==1
    aorta_2 = rois['aorta_2_accept']==1
    liver_1 = rois['liver_1_accept']==1
    liver_2 = rois['liver_2_accept']==1
    t0 = rois['time_1'][0]
    xdata = (
        rois['time_1'][aorta_1] - t0, 
        rois['time_2'][aorta_2] - t0, 
        rois['time_1'][liver_1] - t0,
        rois['time_2'][liver_2] - t0,
    )
    ydata = (
        rois['aorta_1'][aorta_1], 
        rois['aorta_2'][aorta_2], 
        rois['liver_1'][liver_1],
        rois['liver_2'][liver_2],
    )

    # Fit model to data
    model = dc.AortaLiver2scan(
        weight=pars['weight'],
        agent='gadoxetate',
        dose=pars['dose_1'],
        dose2=pars['dose_2'],
        rate=1,
        field_strength=3.0,
        t0=pars['t0'],
        TR=pars['TR'], 
        FA=pars['FA_1'],
        FA2=pars['FA_2'],
        TS=rois['time_1'][1]-t0,
        R10a=1/pars['T1_aorta_1'],
        R10l=1/pars['T1_liver_1'],
        R102a=1/pars['T1_aorta_3'],
        R102l=1/pars['T1_liver_3'],
        H=0.45,
        vol=pars['liver_volume'],
    )
    loss0 = model.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss0)
    model.train(xdata, ydata, xtol=1e-3, verbose=2)
    loss1 = model.cost(xdata, ydata)
    print('Goodness of fit (improvement, %): ', 100*(loss0-loss1)/loss0)

    # Export data
    figure(model, xdata, ydata, os.path.join(path, 'Figs'), name, 
           pars, t0)
    pars = parameters(model, xdata, ydata, t0, pars)
    tools.to_csv(os.path.join(path, 'Pars'), name, pars)


def parameters(model:dc.AortaLiver2scan, xdata, ydata, t0, params):

    time = (xdata[0], xdata[1])
    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[2], ydata[2]

    model.dose2 = 0

    # Compute AUC over 3hrs
    model.tmax = model.BAT+180*60
    t, cb, Cl = model.conc()
    t, R1b, R1l = model.relax()
    AUC_Cb = np.trapezoid(cb, model.t) 
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
    AUC35_Cb = np.trapezoid(cb, model.t)
    AUC35_Cl = np.trapezoid(Cl, t) 

    # Export parameters
    pars = model.export_params()
    pars['AUC_Cb']=['AUC for Cb (0-inf)', 1000*AUC_Cb, 'mM*sec',0]
    pars['AUC_Cl']=['AUC for Cl (0-inf)', 1000*AUC_Cl, 'mM*sec',0]  
    pars['AUC35_Cb']=['AUC for Cb (0-35min)', 1000*AUC35_Cb, 'mM*sec',0]
    pars['AUC35_Cl']=['AUC for Cl (0-35min)', 1000*AUC35_Cl, 'mM*sec',0]  
    pars['RE_R1b']=['RE for R1b at 20min', 100*RE_R1b, '%',0]
    pars['RE_R1l']=['RE for R1l at 20min', 100*RE_R1l, '%',0]
    pars['RE_Sb']=['RE for Sb at 20min', 100*RE_Sb, '%',0]
    pars['RE_Sl']=['RE for Sl at 20min', 100*RE_Sl, '%',0]    

    # MOLLI values     
    pars['T1_1']=['Liver T1-MOLLI at baseline', params['T1_liver_1'], 'sec', 0]
    pars['T1_2']=['Liver T1-MOLLI at 45min', params['T1_liver_2'], 'sec', 0]
    pars['T1_3']=['Liver T1-MOLLI at scan 2', params['T1_liver_3'], 'sec', 0]

    # Timings needed for plotting etc
    pars['t1_MOLLI']=['Time of T1-MOLLI at baseline', params['T1_time_1']/(60*60), 'hrs', 0]
    pars['t2_MOLLI']=['Time of T1-MOLLI at 45min', params['T1_time_2']/(60*60), 'hrs', 0]
    pars['t3_MOLLI']=['Time of T1-MOLLI at scan 2', params['T1_time_3']/(60*60), 'hrs', 0]

    pars['t0']=["Start time first acquisition", (t0+time[0][0])/(60*60), 
                'hrs',0]
    pars['t1']=["End time first acquisition", (t0+time[0][-1])/(60*60), 
                'hrs',0]
    pars['t2']=["Start time second acquisition", (t0+time[1][0])/(60*60), 
                'hrs',0]
    pars['t3']=["End time second acquisition", (t0+time[1][-1])/(60*60), 
                'hrs',0]
    pars['dt1']=["Time step first acquisition", model.TS, 'sec',0]
    pars['dt2']=["Time step second acquisition", model.TS, 'sec',0]

    return pars  


def figure(model:dc.AortaLiver2scan, xdata, ydata, path, name, 
           params, t0):

    t = [
        0, 
        params['T1_time_2']-t0,
        params['T1_time_3']-t0,
    ]
    R1a = [
        1/params['T1_aorta_1'], 
        1/params['T1_aorta_2'],
        1/params['T1_aorta_3'],
    ]
    R1l = [
        1/params['T1_liver_1'], 
        1/params['T1_liver_2'],
        1/params['T1_liver_3'],
    ]

    if not os.path.exists(path):
        os.makedirs(path)

    ya = [dc.signal_ss(model.S0a, R1a[0], model.TR, model.FA),
          dc.signal_ss(model.S0a, R1a[1], model.TR, model.FA),
          dc.signal_ss(model.S02a, R1a[2], model.TR, model.FA)]
    yl = [dc.signal_ss(model.S0l, R1l[0], model.TR, model.FA),
          dc.signal_ss(model.S0l, R1l[1], model.TR, model.FA),
          dc.signal_ss(model.S02l, R1l[2], model.TR, model.FA)]
    test = ((t,ya),(t,yl))
    file = os.path.join(path, name)
    model.plot(xdata, ydata, fname=file + '.png', ref=test, show=False)
    BAT = model.BAT
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+1200], 
               fname=file + '_scan1_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_scan1_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_scan1_win3.png', ref=test, show=False)
    BAT = model.BAT2
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+1200], 
               fname=file + '_scan2_win1.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+600], 
               fname=file + '_scan2_win2.png', ref=test, show=False)
    model.plot(xdata, ydata, xlim=[BAT-20, BAT+160], 
               fname=file + '_scan2_win3.png', ref=test, show=False)

