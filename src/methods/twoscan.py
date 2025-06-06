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


def _data(rois):
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
    return xdata, ydata


def subject_model(data, subj, visit, verbose=0):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    # Define default model
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

    # Personalise model
    xdata, ydata = _data(rois)
    model.train(xdata, ydata, xtol=1e-3, verbose=verbose)
    return model


def save_plots(model, data, subj, visit, path):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    name = subj + '_' + visit
    xdata, ydata = _data(rois)
    t0 = rois['time_1'][0]
    path = os.path.join(path, 'Plots')
    file = os.path.join(path, name)

    t = [
        0, 
        pars['T1_time_2']-t0,
        pars['T1_time_3']-t0,
    ]
    R1a = [
        1/pars['T1_aorta_1'], 
        1/pars['T1_aorta_2'],
        1/pars['T1_aorta_3'],
    ]
    R1l = [
        1/pars['T1_liver_1'], 
        1/pars['T1_liver_2'],
        1/pars['T1_liver_3'],
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


def save_results(model, data, subj, visit, path):

    rois = data['rois'][subj][visit]
    pars = data['pars'][subj][visit]

    xdata, ydata = _data(rois)
    t0 = rois['time_1'][0]
    tb, Sb, tl, Sl = xdata[0], ydata[0], xdata[2], ydata[2]

    time = (xdata[0], xdata[1])
    model.dose2 = 0

    params = tools.export_params(model, tb, Sb, tl, Sl, pars)
    params['T1_3']=['Liver T1-MOLLI at scan 2', pars['T1_liver_3'], 'sec', 0]
    params['t3_MOLLI']=['Time of T1-MOLLI at scan 2', pars['T1_time_3']/(60*60), 'hrs', 0]

    params['t0']=["Start time first acquisition", (t0+time[0][0])/(60*60), 
                'hrs',0]
    params['t1']=["End time first acquisition", (t0+time[0][-1])/(60*60), 
                'hrs',0]
    params['t2']=["Start time second acquisition", (t0+time[1][0])/(60*60), 
                'hrs',0]
    params['t3']=["End time second acquisition", (t0+time[1][-1])/(60*60), 
                'hrs',0]
    params['dt1']=["Time step first acquisition", model.TS, 'sec',0]
    params['dt2']=["Time step second acquisition", model.TS, 'sec',0]

    params = tools.to_tristan_units(params)
    dmrpath = os.path.join(path, 'Results')
    return tools.to_dmr(dmrpath, subj, visit, params)


