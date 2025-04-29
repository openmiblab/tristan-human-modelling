import os
import time

import dcmri as dc
import pydmr

from tristan import tools


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


def compute_vart(datafile, resultspath, 
                 acq_times = [5,10,15,20,25,30,35,40]):

    start = time.time()

    if not os.path.exists(resultspath):
        os.makedirs(resultspath)

    results = []
    data = pydmr.read(datafile, format='nest')
    for subj in data['rois'].keys():
        for visit in data['rois'][subj].keys():
            for tacq in acq_times:
                name = subj + '_' + visit + '_' + str(tacq).zfill(2)
                file = fit_subj(data, subj, visit, resultspath, name, tacq)
                results.append(file)
    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)
    
    print('Calculation time (mins): ', (time.time()-start)/60)


def fit_subj(data, subj, visit, path, name, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]
    
    xdata = (
        rois['time_1'][rois['aorta_1_accept']] - rois['time_1'][0], 
        rois['time_1'][rois['liver_1_accept']] - rois['time_1'][0],
    )
    ydata = (
        rois['aorta_1'][rois['aorta_1_accept']], 
        rois['liver_1'][rois['liver_1_accept']],
    )

    # Truncate data if requested
    if tacq is not None:
        idx0, idx1 = xdata[0]<tacq*60, xdata[1]<tacq*60
        xdata = (xdata[0][idx0], xdata[1][idx1])
        ydata = (ydata[0][idx0], ydata[1][idx1])
    
    # Fit model to data
    model = dc.AortaLiver(

        # Injection parameters
        weight=pars['weight'],
        agent='gadoxetate',
        dose=pars['dose_1'],
        rate=1,

        # Acquisition parameters
        field_strength=3.0,
        t0=pars['t0'],
        TR=pars['TR'], 
        FA=pars['FA_1'],
        TS=rois['time_1'][1]-rois['time_1'][0],

        # Signal parameters
        R10a=1/pars['T1_aorta_1'],
        R10l=1/pars['T1_liver_1'],

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
    pars = parameters(model, xdata, ydata, pars)
    study = visit if tacq is None else visit + '_' + str(tacq).zfill(2)
    dmrpath = os.path.join(path, 'Results')
    return tools.to_dmr(dmrpath, subj, study, name, pars)


def parameters(model:dc.AortaLiver, xdata, ydata, params):

    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[1], ydata[1]
    pars = tools.export_params(model, tb, Sb, tl, Sl, params) 
    return tools.to_tristan_units(pars)


def figure(model:dc.AortaLiver, xdata, ydata, path, name, params, t0):
    
    t = [
        0, 
        params['T1_time_2']-t0,
    ]
    R1a = [
        1/params['T1_aorta_1'], 
        1/params['T1_aorta_2'],
    ]
    R1l = [
        1/params['T1_liver_1'], 
        1/params['T1_liver_2'],
    ]

    if not os.path.exists(path):
        os.makedirs(path)
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

