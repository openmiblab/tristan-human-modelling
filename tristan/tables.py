import os
import pandas as pd
import numpy as np

from tristan import calc





def cases(results, folder):

    path = os.path.join(results, folder, 'Tables')
    if not os.path.exists(path):
        os.makedirs(path)
        
    # Get data
    df = pd.read_csv(
        os.path.join(results, folder, 'Analysis', 'parameters_rep.csv')
    )
    src = os.path.join(results, folder)
    df['group'] = calc.lookup(src, df.parameter.values, 'group')
    df['description'] = calc.lookup(src, df.parameter.values, 'description')
    df['unit'] = calc.lookup(src, df.parameter.values, 'unit') 

    # Get effect sizes
    ef = pd.read_csv(
        os.path.join(results, folder, 'Analysis', 'effect_size.csv')
    )
    ef['group'] = calc.lookup(src, ef.parameter.values, 'group')
    ef['description'] = calc.lookup(src, ef.parameter.values, 'description')
    ef['visit'] = 'change (%)'
    ef['unit'] = calc.lookup(src, ef.parameter.values, 'unit')
    ef['stdev'] = np.nan

    # Concatenate
    ef = ef[['subject','visit','parameter','description','value','unit','stdev','group']]
    df = df[['subject','visit','parameter','description','value','unit','stdev','group']]
    df = pd.concat([df, ef])
    dfa = df[df.group=='MRI - aorta']
    df = df[df.group=='MRI - liver']

    for i, subject in enumerate(df.subject.unique()):

        subj = str(subject).zfill(3)
        dfs = df[df.subject==subject]
        dfas = dfa[dfa.subject==subject]
        
        # Liver values
        pivot = pd.pivot_table(dfs, values='value', columns='visit', 
                                index=['description','unit'])
        cols = pivot.columns.tolist()
        if len(cols)>1:
            pivot = pivot[['control', 'drug', 'change (%)']]
        header = ['Biomarker', 'Units'] + list(pivot.columns)
        table = []
        for row in pivot.index:
            rvals = list(np.around(pivot.loc[row,:].values,2))
            table.append([row[0], row[1]] + rvals)
        file = os.path.join(path, subj+'_liver.csv')
        pd.DataFrame(table, columns=header).to_csv(file, index=False)

        # Aorta values
        pivot = pd.pivot_table(dfas, values='value', columns='visit', 
                                index=['description','unit'])
        cols = pivot.columns.tolist()
        if len(cols)>1:
            pivot = pivot[['control', 'drug', 'change (%)']]
        header = ['Biomarker', 'Units'] + list(pivot.columns)
        table = []
        for row in pivot.index:
            rvals = list(np.around(pivot.loc[row,:].values,2))
            table.append([row[0], row[1]] + rvals)
        file = os.path.join(path, subj+'_aorta.csv')
        pd.DataFrame(table, columns=header).to_csv(file, index=False)


def reference(results, folder):

    path = os.path.join(results, folder, 'Tables')
    if not os.path.exists(path):#
        os.makedirs(path)

    for visit in ['control', 'drug']:
        for group in ['aorta', 'liver']:
            file = os.path.join(results, folder, 'Analysis', 'difference_with_reference.csv')
            df = pd.read_csv(file).set_index('Biomarker')
            df = df[df.visit==visit]
            df = df[df.group=='MRI - '+group]
            df.drop(columns=['visit','group'], inplace=True)
            header = [df.index.name] + list(df.columns)
            table = []
            for row in df.index:
                table.append([row] + list(df.loc[row,:]))
            file = os.path.join(path, f'reference_{group}_{visit}.csv')
            pd.DataFrame(table, columns=header).to_csv(file, index=False)


def pairwise_diff(results, folder):

    path = os.path.join(results, folder, 'Tables')
    if not os.path.exists(path):
        os.makedirs(path)

    data = pd.read_csv(os.path.join(results, folder, 'Analysis', 'avr_95CI.csv'))
    b_avr = data['mean control'].values
    r_avr = data['mean drug'].values
    c_avr = data['mean effect'].values
    b_err = data['95%CI control'].values
    r_err = data['95%CI drug'].values
    c_err = data['95%CI effect'].values

    # Update output array
    data['control'] = [str(b_avr[i]) + ' (' + str(b_err[i]) + ') ' 
                         for i in range(b_avr.size)]
    data['drug'] = [str(r_avr[i]) + ' (' + str(r_err[i]) + ') ' 
                         for i in range(r_avr.size)]
    data['change (%)'] = [str(c_avr[i]) + ' (' + str(c_err[i]) + ') ' 
                            for i in range(c_avr.size)]
    
    data = data[['description','group','unit', 'control', 'drug', 'change (%)']]
    data = data.rename(columns={'description': "Biomarker", "unit": "Units"})

    liver = data[data.group=='MRI - liver']
    liver = liver.drop(columns='group')
    liver = liver.sort_values('Biomarker', ascending=True)
    liver.to_csv(os.path.join(path, 'liver_pairwise.csv'), index=False)

    aorta = data[data.group=='MRI - aorta']
    aorta = aorta.drop(columns='group')
    aorta = aorta.sort_values('Biomarker', ascending=True)
    aorta.to_csv(os.path.join(path, 'aorta_pairwise.csv'), index=False)



def pairwise_stats(results, folder):

    path = os.path.join(results, folder, 'Tables')
    if not os.path.exists(path):
        os.makedirs(path)

    output = pd.read_csv(os.path.join(results, folder, 'Analysis', '_output_ttest.csv'))

    output.drop(columns=['Contrast', 'A', 'B', 'Paired', 'Parametric', 'T', 
                         'dof', 'alternative'], inplace=True)
    output = output[['group','description', 'p-unc', 'BF10', 'odds-ratio']]
    output.rename(columns={'description': "Biomarker", 
                           "p-unc": "p-value", 'BF10': 'Bayes Factor', 
                           'odds-ratio': 'Odds Ratio'}, inplace=True)
    
    output = output.astype({'Bayes Factor': 'float32'})
    output.loc[:,'p-value'] = np.around(output['p-value'].values, 5)
    output.loc[:,'Bayes Factor'] = np.around(output['Bayes Factor'].values, 2)
    output.loc[:,'Odds Ratio'] = np.around(output['Odds Ratio'].values, 2)

    output = output[['Biomarker', 'group', "p-value", 'Bayes Factor', 'Odds Ratio']]
    
    liver = output[output.group=='MRI - liver']
    liver.drop(columns='group', inplace=True)
    liver = liver.sort_values('Biomarker', ascending=True)
    liver.to_csv(os.path.join(path, 'liver_ttest.csv'), index=False)

    aorta = output[output.group=='MRI - aorta']
    aorta.drop(columns='group', inplace=True)
    aorta = aorta.sort_values('Biomarker', ascending=True)
    aorta.to_csv(os.path.join(path, 'aorta_ttest.csv'), index=False)


# def ttest(resultspath, folder):

#     path = os.path.join(resultspath, folder, 'Tables')
#     if not os.path.exists(path):#
#         os.makedirs(path)

#     df = pd.read_csv(os.path.join(resultspath, folder, 'Analysis', 'pairwise_statistics.csv')).set_index('Biomarker')
#     df = df[df.group=='MRI - liver']
#     df.drop(columns='group', inplace=True)
#     df.to_csv(os.path.join(path, 'liver_ttest.csv'))

#     df = pd.read_csv(os.path.join(resultspath, folder, 'Analysis', 'pairwise_statistics.csv')).set_index('Biomarker')
#     df = df[df.group=='MRI - aorta']
#     df.drop(columns='group', inplace=True)
#     df.to_csv(os.path.join(path, 'aorta_ttest.csv'))

    # df = pd.read_csv(os.path.join(resultspath, folder, 'Analysis', 'pairwise_differences.csv')).set_index('Biomarker')
    # df = df[df.group=='MRI - liver']
    # df.drop(columns='group', inplace=True)
    # df.to_csv(os.path.join(path, 'liver_pairwise.csv'))

    # df = pd.read_csv(os.path.join(resultspath, folder, 'Analysis', 'pairwise_differences.csv')).set_index('Biomarker')
    # df = df[df.group=='MRI - aorta']
    # df.drop(columns='group', inplace=True)
    # df.to_csv(os.path.join(path, 'aorta_pairwise.csv'))


