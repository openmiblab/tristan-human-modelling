import os
import pickle
import pandas as pd

# https://bioxydynltd.sharepoint.com/:f:/s/WP2-TRISTAN/Enz3hJ8UPddEgZxXdaspo5YB7GkeLKsJB5k1Xd2mmr8jPw?e=5%3avvdrgU&at=9


def read(subj):

    const = pd.read_excel(subj, sheet_name='const')
    const.set_index('name', inplace=True)
    dyn1 = pd.read_excel(subj, sheet_name='dyn1')
    dyn2 = pd.read_excel(subj, sheet_name='dyn2')
    molli1 = pd.read_excel(subj, sheet_name='MOLLI1')
    molli2 = pd.read_excel(subj, sheet_name='MOLLI2')
    molli3 = pd.read_excel(subj, sheet_name='MOLLI3')
    dyn1.sort_values('time', inplace=True)
    dyn2.sort_values('time', inplace=True)
    molli1.sort_values('time', inplace=True)
    molli2.sort_values('time', inplace=True)
    molli3.sort_values('time', inplace=True)
    if 'aorta_valid' not in dyn1:
        dyn1['aorta_valid'] = 1
    if 'liver_valid' not in dyn1:
        dyn1['liver_valid'] = 1
    if 'kidney_valid' not in dyn1:
        dyn1['kidney_valid'] = 1
    if 'portal_valid' not in dyn1:
        dyn1['portal_valid'] = 1
    if 'aorta_valid' not in dyn2:
        dyn2['aorta_valid'] = 1
    if 'liver_valid' not in dyn2:
        dyn2['liver_valid'] = 1
    if 'kidney_valid' not in dyn2:
        dyn2['kidney_valid'] = 1
    if 'portal_valid' not in dyn2:
        dyn2['portal_valid'] = 1
    t0 = dyn1.time.values[0]
    
    data = {
        'time1': dyn1.time.values-t0, 
        'fa1': dyn1.fa.values, 
        'aorta1': dyn1.aorta.values, 
        'liver1': dyn1.liver.values, 
        'kidney1': dyn1.kidney.values,
        'portal1': dyn1['portal-vein'].values, 
        'aorta_valid1': dyn1.aorta_valid.values, 
        'liver_valid1': dyn1.liver_valid.values, 
        'portal_valid1': dyn1.portal_valid.values,
        'kidney_valid1': dyn1.kidney_valid.values,
        'time2': dyn2.time.values-t0, 
        'fa2': dyn2.fa.values, 
        'aorta2': dyn2.aorta.values, 
        'liver2': dyn2.liver.values, 
        'kidney2': dyn2.kidney.values,
        'portal2': dyn2['portal-vein'].values, 
        'aorta_valid2': dyn2.aorta_valid.values, 
        'liver_valid2': dyn2.liver_valid.values, 
        'portal_valid2': dyn2.portal_valid.values,
        'kidney_valid2': dyn2.kidney_valid.values,
        'T1time1': molli1.time.values[0]-t0, 
        'T1aorta1': molli1.aorta.values[0], 
        'T1liver1': molli1.liver.values[0], 
        'T1kidney1': molli1['kidney-parenchyma'].values[0], 
        'T1portal1': molli1['portal vein'].values[0],
        'T1time2': molli2.time.values[0]-t0, 
        'T1aorta2': molli2.aorta.values[0], 
        'T1liver2': molli2.liver.values[0], 
        'T1kidney2': molli2['kidney-parenchyma'].values[0], 
        'T1portal2': molli2['portal vein'].values[0],
        'T1time3': molli3.time.values[0]-t0, 
        'T1aorta3': molli3.aorta.values[0], 
        'T1liver3': molli3.liver.values[0], 
        'T1kidney3': molli3['kidney-parenchyma'].values[0], 
        'T1portal3': molli3['portal vein'].values[0],
        'weight': const.at['weight', 'value'], 
        'dose1': const.at['dose1', 'value'], 
        'dose2': const.at['dose2', 'value'], 
        'baseline': const.at['baseline', 'value'], 
        'liver_volume': const.at['liver-volume-mm3','value']/1000,
        'kidney_volume': const.at['kidney-volume-mm3','value']/1000,
        't0': t0,
    }

    # Format data
    xa = (data['aorta_valid1']==1, 
          data['aorta_valid2']==1)
    xl = (data['liver_valid1']==1, 
          data['liver_valid2']==1)
    data['xdata'] = (
             data['time1'][xa[0]], 
             data['time2'][xa[1]],
             data['time1'][xl[0]], 
             data['time2'][xl[1]])
    data['ydata'] = (
             data['aorta1'][xa[0]], 
             data['aorta2'][xa[1]], 
             data['liver1'][xl[0]], 
             data['liver2'][xl[1]])
    data['tR1'] = [0, data['T1time2'], data['T1time3']]
    data['R1l'] = [
           1000.0/data['T1liver1'], 
           1000.0/data['T1liver2'],
           1000.0/data['T1liver3']]
    data['R1a'] = [
           1000.0/data['T1aorta1'], 
           1000.0/data['T1aorta2'],
           1000.0/data['T1aorta3']]

    return data



datapath = os.path.dirname(os.path.abspath(__file__))

data_list = []
for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
    visitdatapath = os.path.join(datapath, visit)
    for s in os.listdir(visitdatapath):
        subj = os.path.join(visitdatapath, s)
        subj_data = read(subj)
        vals = {
            'visit': visit,
            'subject': s[:3],
            'time1aorta': subj_data['xdata'][0],
            'time2aorta': subj_data['xdata'][1],
            'time1liver': subj_data['xdata'][2],
            'time2liver': subj_data['xdata'][3],
            'signal1aorta': subj_data['ydata'][0],
            'signal2aorta': subj_data['ydata'][1],
            'signal1liver': subj_data['ydata'][2],
            'signal2liver': subj_data['ydata'][3],
            'tR1':  subj_data['tR1'],
            'R1a':  subj_data['R1a'],
            'R1l':  subj_data['R1l'],
            'tacq':  subj_data['time1'][2]-subj_data['time1'][1],
            'tacq2': subj_data['time2'][2]-subj_data['time2'][1],
            'weight':  subj_data['weight'],
            'agent': 'gadoxetate',
            'dose':  [subj_data['dose1'], subj_data['dose2']],
            'rate':  1,
            'field_strength':  3.0,
            't0':  subj_data['baseline'],
            'TR':  3.71/1000.0,
            'FA':  15,
            'R10b': subj_data['R1a'][0],
            'R10l': subj_data['R1l'][0],
            'R102b': subj_data['R1a'][2],
            'R102l': subj_data['R1l'][2],
            'Hct':  0.45,
            'vol':  subj_data['liver_volume'],
        }
        data_list.append(vals)

with open(os.path.join(datapath, 'tristan_rifampicin.pkl'), 'wb') as fp:
    pickle.dump(data_list, fp, protocol=pickle.HIGHEST_PROTOCOL)





