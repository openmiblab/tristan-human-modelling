import os
import time

import dcmri as dc
import pydmr

from methods import tools



def compute(datafile, resultspath):

    start = time.time()

    if not os.path.exists(resultspath):
        os.makedirs(resultspath)

    results = []
    data = pydmr.read(datafile, format='nest')
    for subj in data['rois'].keys():
        for visit in data['rois'][subj].keys():
            name = subj + '_' + visit
            file = fit_subj(data, subj, visit, resultspath, name)
            results.append(file)
    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)

    print('Calculation time (mins): ', (time.time()-start)/60)


def fit_subj(data, subj, visit, path, name):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    xdata = (
        rois['time_1'][rois['aorta_1_accept']] - rois['time_1'][0], 
        rois['time_2'][rois['aorta_2_accept']] - rois['time_1'][0], 
        rois['time_1'][rois['liver_1_accept']] - rois['time_1'][0],
        rois['time_2'][rois['liver_2_accept']] - rois['time_1'][0],
    )
    ydata = (
        rois['aorta_1'][rois['aorta_1_accept']], 
        rois['aorta_2'][rois['aorta_2_accept']], 
        rois['liver_1'][rois['liver_1_accept']],
        rois['liver_2'][rois['liver_2_accept']],
    )

    # Fit model to data
    model = dc.AortaLiver2scan(

        # Injection parameters
        weight=pars['weight'],
        agent='gadoxetate',
        dose=pars['dose_1'],
        dose2=pars['dose_2'],
        rate=1,

        # Acquisition parameters
        field_strength=3.0,
        t0=pars['t0'],
        TR=pars['TR'], 
        FA=pars['FA_1'],
        FA2=pars['FA_2'],
        TS=rois['time_1'][1]-rois['time_1'][0],

        # Signal parameters
        R10a=1/pars['T1_aorta_1'],
        R10l=1/pars['T1_liver_1'],
        R102a=1/pars['T1_aorta_3'],
        R102l=1/pars['T1_liver_3'],

        # Tissue parameters
        vol=pars['liver_volume'],
    )
    loss0 = model.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss0)
    model.train(xdata, ydata, xtol=1e-3, verbose=2)
    loss1 = model.cost(xdata, ydata)
    print('Goodness of fit (improvement, %): ', 100*(loss0-loss1)/loss0)

    # Export data
    pngpath =  os.path.join(path, 'Plots')
    figure(model, xdata, ydata, pngpath, name, pars, rois['time_1'][0])
    pars = parameters(model, xdata, ydata, rois['time_1'][0], pars)
    dmrpath = os.path.join(path, 'Results')
    return tools.to_dmr(dmrpath, subj, visit, name, pars)


def parameters(model:dc.AortaLiver2scan, xdata, ydata, t0, params):

    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[2], ydata[2]

    time = (xdata[0], xdata[1])
    model.dose2 = 0

    pars = tools.export_params(model, tb, Sb, tl, Sl, params)
    pars['T1_3']=['Liver T1-MOLLI at scan 2', params['T1_liver_3'], 'sec', 0]
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

    return tools.to_tristan_units(pars)


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

