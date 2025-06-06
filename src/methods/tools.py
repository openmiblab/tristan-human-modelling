import os

import numpy as np

import pydmr


AORTA_PARS = ['RE_Sb', 'RE_R1b', 'S02a', 'BAT','CO','Thl','Dhl',
              'To','Eb','Eo','Toe','BAT2', 
              'AUC_R1b','AUC_Cb', 'AUC35_R1b','AUC35_Cb']

LABEL = {
    'AUC_Cb': 'AUC(b)',
    'AUC_Cl': 'AUC(l)',
    'AUC_R1b': 'AUR1(b)',
    'AUC_R1l': 'AUR1(l)',
    'AUC35_Cb': 'AUC(30,b)',
    'AUC35_Cl': 'AUC(30,l)',
    'AUC35_R1b': 'AUR1(30,b)',
    'AUC35_R1l': 'AUR1(30,l)',
    'BAT': 'BAT(1)',
    'BAT2': 'BAT(2)',
    'CL': 'Cl',
    'CO': 'CO',
    'De': 'TTD(e)',
    'Dhl': 'TTD(hl)',
    'dt1': 'dt(1)',
    'dt2': 'dt(2)',
    'Eb': 'E',
    'Eo': 'E(o)',
    'Kbh': 'K(bh)',
    'kbh': 'k(bh)',
    'kbh_f': 'k(bh,f)',
    'kbh_i': 'k(bh,i)',
    'Khe': 'K(he)',
    'khe': 'k(he)',
    'khe_f': 'k(he,f)',
    'khe_i': 'k(he,i)',
    'RE_R1b': 'RER1(b)',
    'RE_R1l': 'RER1(l)',
    'RE_Sb': 'RES(b)',
    'RE_Sl': 'RES(l)',
    'S02a': 'S0(2,a)',
    'S02l': 'S0(2,l)',
    't0': 't0',
    't1': 't1',
    'T1_1': 'T1(1)',
    'T1_2': 'T1(45)',
    'T1_3': 'T1(2)',
    't1_MOLLI': 't1_MOLLI',
    't2': 't2',
    't2_MOLLI': 't2_MOLLI',
    't3': 't3',
    't3_MOLLI': 't3_MOLLI',
    'Te': 'MTT(e)',
    'Th': 'MTT(h)',
    'Th_f': 'MTT(h,f)',
    'Th_i': 'MTT(h,i)',
    'Thl': 'MTT(hl)',
    'To': 'MTT(o)',
    'Toe': 'MTT(o,e)',
    've': 'v(e)',
    'H': 'Hct',
}


def export_params(model, tb, Sb, tl, Sl, params):

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

    # Timings needed for plotting etc
    pars['t1_MOLLI']=['Time of T1-MOLLI at baseline', params['T1_time_1']/(60*60), 'hrs', 0]
    pars['t2_MOLLI']=['Time of T1-MOLLI at 45min', params['T1_time_2']/(60*60), 'hrs', 0]
  
    return pars  



def to_dmr(path, subj, study, pars):

    # Build data dictionary
    dmr = {
        'data': {},
        'pars': {},   
        'sdev': {},
        'columns': ['group', 'label'],
    }
    for key, val in pars.items():
        dmr['data'][key] = [val[0], val[2], 'float']
        dmr['pars'][subj, study, key] = val[1]
        dmr['sdev'][subj, study, key] = val[3]

    # Append group and label to the data dictionary
    for p in dmr['data']:
        if p in AORTA_PARS:
            dmr['data'][p].append('MRI - aorta')
        else:
            dmr['data'][p].append('MRI - liver')
        dmr['data'][p].append(LABEL[p])

    # Save as dmr file
    name = subj + '_' + study
    file = os.path.join(path, name + '.dmr')
    pydmr.write(file, dmr)
    return file



def to_tristan_units(pars):
  
    slow_time = ['BAT','BAT2','Toe','Th','Th_i','Th_f']
    for p in pars:
        if pars[p][2] == '':
            pars[p][1:] = [pars[p][1]*100, '%', pars[p][3]*100]
        if pars[p][2] == 'mL/cm3':
            pars[p][1:] = [pars[p][1]*100, 'mL/100cm3', pars[p][3]*100]
        if pars[p][2] == 'mL/sec/cm3':
            pars[p][1:] = [pars[p][1]*6000, 'mL/min/100cm3', pars[p][3]*6000]
        if p in slow_time:
            pars[p][1:] = [pars[p][1]/60, 'min', pars[p][3]/60]

    if 'CO' in pars:
        pars['CO'][1:] = [pars['CO'][1]*60/1000, 'L/min', pars['CO'][3]*60/1000]
    if 'CL' in pars:
        pars['CL'][1:] = [pars['CL'][1]*60/1000, 'L/min', pars['CL'][3]*60/1000]

    return pars


