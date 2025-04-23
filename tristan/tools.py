import os
import time

import pandas as pd


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
}


    

def build_master(resultspath, vart=False):
    path = os.path.join(resultspath, 'Pars')
    filenames = os.listdir(path)
    output = None
    for filename in filenames:
        subj_file = os.path.join(path, filename)
        pars = pd.read_csv(subj_file)
        pars['subject'] = filename[:3]
        pars['visit'] = 'control' if 'control' in filename else 'drug'
        group = []
        for p in pars.parameter.values:
            if p in AORTA_PARS:
                group.append('MRI - aorta')
            else:
                group.append('MRI - liver')
        pars['group'] = group
        pars['label'] = pars['parameter'].map(LABEL)
        if vart:
            pars['tacq'] = float(filename[-6:-4])
        if output is None:
            output = pars
        else:
            output = pd.concat([output, pars])

    # Export data for analysis (final format)
    labels = ['subject','visit','parameter','value','stdev']
    if vart:
        labels.append('tacq')
    output[labels].to_csv(os.path.join(resultspath, 'output_data.csv'), index=False)

    # Export dictionary for analysis
    pars = pars[['parameter','name','unit','group','label']]
    pars.to_csv(os.path.join(resultspath, 'output_data_dict.csv'),index=False)


# def build_master_vart(resultspath):
#     path = os.path.join(resultspath, 'Pars')
#     filenames = os.listdir(path)
#     output = None
#     for filename in filenames:
#         subj_file = os.path.join(path, filename)
#         pars = pd.read_csv(subj_file)
#         pars['subject'] = filename[:3]
       
#         pars['visit'] = 'control' if 'control' in filename else 'drug'
#         structure = []
#         for p in pars.parameter.values:
#             if p in AORTA_PARS:
#                 structure.append('MRI - aorta')
#             else:
#                 structure.append('MRI - liver')
#         pars['group'] = structure
#         if output is None:
#             output = pars
#         else:
#             output = pd.concat([output, pars])
#         pars['tacq'] = float(filename[-6:-4])

#     # Export data for analysis (final format)
#     output = output[['subject','visit','parameter','value','stdev', 'tacq']]
#     output.to_csv(os.path.join(resultspath, 'output_data.csv'), index=False)

#     # Export dictionary for analysis
#     pars = pars[['parameter','name','unit','group']]
#     pars.to_csv(os.path.join(resultspath, 'output_data_dict.csv'),index=False)



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



def to_csv(path, name, pars):
    pars = to_tristan_units(pars)
    if not os.path.exists(path):
        os.makedirs(path)
    file = os.path.join(path, name + '.csv')
    cols = ['parameter', 'name', 'value', 'unit', 'stdev']
    pars = [[key] + val for key, val in pars.items()]
    df = pd.DataFrame(pars, columns=cols)
    path = os.path.dirname(file)
    if not os.path.isdir(path):
        os.makedirs(path)
    try:
        df.to_csv(file, index=False)
    except:
        print("Can't write to file ", file)
        print("Please close the file before saving data.") 


