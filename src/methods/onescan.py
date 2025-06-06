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

            # Train model
            model = subject_model(data, subj, visit, verbose=2)

            # Save results
            save_plots(model, data, subj, visit, resultspath)
            file = save_results(model, data, subj, visit, resultspath)

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
                
                # Train model
                model = subject_model(data, subj, visit, verbose=2, tacq=tacq)

                # Save results
                save_plots(model, data, subj, visit, resultspath, tacq=tacq)
                file = save_results(model, data, subj, visit, resultspath, tacq=tacq)

                results.append(file)
    file = os.path.join(resultspath, 'all_results')
    pydmr.concat(results, file)
    
    print('Calculation time (mins): ', (time.time()-start)/60)


def _data(rois):

    xdata = (
        rois['time_1'][rois['aorta_1_accept']] - rois['time_1'][0], 
        rois['time_1'][rois['liver_1_accept']] - rois['time_1'][0],
    )
    ydata = (
        rois['aorta_1'][rois['aorta_1_accept']], 
        rois['liver_1'][rois['liver_1_accept']],
    )
    return xdata, ydata


def subject_model(data, subj, visit, verbose=0, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

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

    # Personalise model

    # Truncate data if requested
    xdata, ydata = _data(rois)
    if tacq is not None:
        idx0, idx1 = xdata[0]<tacq*60, xdata[1]<tacq*60
        xdata = (xdata[0][idx0], xdata[1][idx1])
        ydata = (ydata[0][idx0], ydata[1][idx1])

    # TRain
    model.train(xdata, ydata, xtol=1e-3, verbose=verbose)
    return model


def save_plots(model, data, subj, visit, path, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    study = visit if tacq is None else visit + '_' + str(tacq).zfill(2)
    name = subj + '_' + study
    xdata, ydata = _data(rois)
    t0 = rois['time_1'][0]
    path = os.path.join(path, 'Plots')
    file = os.path.join(path, name)
    
    t = [
        0, 
        pars['T1_time_2']-t0,
    ]
    R1a = [
        1/pars['T1_aorta_1'], 
        1/pars['T1_aorta_2'],
    ]
    R1l = [
        1/pars['T1_liver_1'], 
        1/pars['T1_liver_2'],
    ]

    if not os.path.exists(path):
        os.makedirs(path)
    
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


def save_results(model, data, subj, visit, path, tacq=None):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]
    xdata, ydata = _data(rois)

    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[1], ydata[1]
    params = tools.export_params(model, tb, Sb, tl, Sl, pars)

    params = tools.to_tristan_units(params)
    dmrpath = os.path.join(path, 'Results')

    study = visit if tacq is None else visit + '_' + str(tacq).zfill(2)
    return tools.to_dmr(dmrpath, subj, study, params)