import os

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


def build_master(resultspath):
    path = os.path.join(resultspath, 'Results')
    filenames = os.listdir(path)

    # combine dmr files
    output_data = {}
    output_pars = {}
    output_sdev = {}
    for filename in filenames:
        subj_file = os.path.join(path, filename)
        dmr = pydmr.read(subj_file)
        output_data = output_data | dmr['data']
        output_pars = output_pars | dmr['pars']
        output_sdev = output_sdev | dmr['sdev']
    
    # Append group and label to the data dictionary
    for p in output_data:
        if p in AORTA_PARS:
            output_data[p].append('MRI - aorta')
        else:
            output_data[p].append('MRI - liver')
        output_data[p].append(LABEL[p])

    dmr = {
        'data': output_data,
        'pars': output_pars,
        'sdev': output_sdev,
        'columns': ['group', 'label'],
    }
    dmr_file = os.path.join(resultspath, 'all_results')
    pydmr.write(dmr_file, dmr)


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


def to_dmr(path, subj, study, name, pars):
    pars = to_tristan_units(pars)
    data_dict = {}
    dmr_pars = {}
    dmr_sdev = {}
    for key, val in pars.items():
        data_dict[key] = [val[0], val[2], 'float']
        dmr_pars[subj, study, key] = val[1]
        dmr_sdev[subj, study, key] = val[3]
    dmr = {
        'data': data_dict,
        'pars': dmr_pars,   
        'sdev': dmr_sdev,
    }
    file = os.path.join(path, name + '.dmr')
    pydmr.write(file, dmr)


