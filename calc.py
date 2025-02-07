import os
import math
import pandas as pd
import numpy as np
import pingouin as pg


def _derive_effect_sizes(output):
    for subject in output.subject.unique():
        df_subj = output[output.subject==subject]
        for structure in df_subj.structure.unique():
            df_struct = df_subj[df_subj.structure==structure]
            for par in df_struct.parameter.unique():
                df_par = df_struct[df_struct.parameter==par]
                if len(df_par)==2:
                    visits = df_par.visit.unique()
                    v0 = df_par[df_par.visit==visits[0]].value.values[0]
                    v1 = df_par[df_par.visit==visits[1]].value.values[0]
                    effect = 100*np.divide(v1-v0, v0)
                    name = df_par[df_par.visit==visits[0]].name.values[0]
                    unit = df_par[df_par.visit==visits[0]].unit.values[0]
                    output.loc[len(output)] = [
                        subject, 'change (%)', structure, name, effect, 
                        unit, 0, par + ' effect size']
    return output


def derive_pars(src, ref=False):
    output_file = os.path.join(src, 'parameters.pkl')
    output = pd.read_pickle(output_file)
    todrop = ['BAT','BAT1','BAT2','S02a','S02']
    todrop += [
        'S02l','t0','t1','t2','t3','dt1','dt2',
        'khe_min','khe_max', 'khe_var',
        'kbh_min','kbh_max', 'kbh_var',
        'Khe_i','Khe_f','Kbh_i','Kbh_f', 'Th_i', 'Th_f']
    todrop += [v + ' effect size' for v in todrop]
    df = output[output.parameter.isin(todrop)]
    output.drop(index=df.index, inplace=True)
    output.reset_index(drop=True, inplace=True)
    if ref:
        output.to_csv(os.path.join(src, 'reference.csv'))
        output.to_pickle(os.path.join(src, 'reference.pkl'))
    else:
        output = _derive_effect_sizes(output)
        output.to_csv(os.path.join(src, 'parameters_rep.csv'))
        output.to_pickle(os.path.join(src, 'parameters_rep.pkl'))



def desc_stats(src):
    #output_file = os.path.join(src, 'parameters_ext.pkl')
    output_file = os.path.join(src, 'parameters_rep.pkl')
    output = pd.read_pickle(output_file)
    visits = output.visit[output.visit!='change (%)'].unique()
    v0 = output[output.visit==visits[0]]
    v0 = pd.pivot_table(v0, values='value', index='subject', 
                        columns='parameter')
    v1 = output[output.visit==visits[1]]
    v1 = pd.pivot_table(v1, values='value', index='subject', 
                        columns='parameter')
    ef = output[output.visit=='change (%)']
    ef = pd.pivot_table(ef, values='value', index='subject', 
                        columns='parameter')
    v0.to_csv(os.path.join(src, '_table_'+visits[0]+'.csv'))
    v1.to_csv(os.path.join(src, '_table_'+visits[1]+'.csv'))
    ef.to_csv(os.path.join(src, '_table_effect.csv'))
    # Calculate stats
    bstats = v0.describe()
    bstats = bstats[['khe', 'kbh']].round(2)
    bstats = bstats.rename(columns={
        "khe": "khe "+visits[0]+" (mL/min/100mL)", 
        "kbh": "kbh "+visits[0]+" (mL/min/100mL)"})
    rstats = v1.describe()
    rstats = rstats[['khe', 'kbh']].round(2)
    rstats = rstats.rename(columns={
        "khe": "khe "+visits[1]+" (mL/min/100mL)", 
        "kbh": "kbh "+visits[1]+" (mL/min/100mL)"})
    estats = ef.describe()
    estats = estats[['khe effect size', 'kbh effect size']].round(1)
    estats = estats.rename(columns={
        "khe effect size": "khe effect size (%)", 
        "kbh effect size": "kbh effect size (%)"})
    stats = pd.concat([estats.T, bstats.T, rstats.T])
    stats=stats.reindex([
        "khe effect size (%)",
        "khe "+visits[0]+" (mL/min/100mL)",
        "khe "+visits[1]+" (mL/min/100mL)",
        "kbh effect size (%)",
        "kbh "+visits[0]+" (mL/min/100mL)",
        "kbh "+visits[1]+" (mL/min/100mL)",
        ])
    stats.to_csv(os.path.join(src, '_table_k_stats.csv'))
    stats.to_pickle(os.path.join(src, '_table_k_stats.pkl'))


def ttest(src):
    
    df = pd.read_pickle(os.path.join(src, 'parameters_rep.pkl'))
    df.fillna({'unit': 'None'}, inplace=True)
    visits = df.visit[df.visit!='change (%)'].unique()

    # Get mean, sdev and count of all variables
    avr = pd.pivot_table(df, values='value', columns='visit', 
                         index=['structure','name','unit'], aggfunc='mean')
    std = pd.pivot_table(df, values='value', columns='visit', 
                         index=['structure','name','unit'], aggfunc='std')
    cnt = pd.pivot_table(df, values='value', columns='visit', 
                         index=['structure','name','unit'], aggfunc='count')

    avr = avr.sort_values(['structure','name']) 
    std = std.sort_values(['structure','name'])
    cnt = cnt.sort_values(['structure','name'])

    # Calculate 95% CI intervals
    b_avr = around_sig(avr[visits[0]].values, 3)
    r_avr = around_sig(avr[visits[1]].values, 3)
    c_avr = around_sig(avr['change (%)'].values, 3)
    b_err = around_sig(
        1.96*std[visits[0]].values / np.sqrt(cnt[visits[0]].values), 2)
    r_err = around_sig(
        1.96*std[visits[1]].values / np.sqrt(cnt[visits[1]].values), 2)
    c_err = around_sig(
        1.96*std['change (%)'].values / np.sqrt(cnt['change (%)'].values), 2)

    # Perform t-tests and save if df output
    output = None
    df = df[df.visit != 'change (%)']
    for struct in df.structure.unique():
        dfs = df[df.structure==struct]
        for par in dfs.name.unique():
            dfp = dfs[dfs.name==par]
            stats = pg.pairwise_tests(data=dfp, 
                    dv='value', within='visit', subject='subject', 
                    return_desc=False, effsize='odds-ratio')
            stats['structure'] = [struct]
            stats['name'] = [par]
            if output is None:
                output = stats
            else:
                output = pd.concat([output, stats])
    output = output.sort_values(['structure','name'])

    # Update output array
    output['unit'] = [i[2] for i in avr.index.tolist()]
    output[visits[0]] = [str(b_avr[i]) + ' (' + str(b_err[i]) + ') ' 
                         for i in range(b_avr.size)]
    output[visits[1]] = [str(r_avr[i]) + ' (' + str(r_err[i]) + ') ' 
                         for i in range(r_avr.size)]
    output['change (%)'] = [str(c_avr[i]) + ' (' + str(c_err[i]) + ') ' 
                            for i in range(c_avr.size)]
    output.reset_index(drop=True, inplace=True)
    output.drop(columns=['Contrast', 'A', 'B', 'Paired', 'Parametric', 'T', 
                         'dof', 'alternative'], inplace=True)
    output = output[['structure','name', 'unit', visits[0], visits[1], 
                     'change (%)', 'p-unc','BF10', 'odds-ratio']]
    output.rename(columns={"name": "Biomarker", "unit": "Units", 
                           "p-unc": "p-value", 'BF10': 'Bayes Factor', 
                           'odds-ratio': 'Odds Ratio'}, inplace=True)
    output = output.sort_values(['structure','p-value'], ascending=True)
    output.set_index('Biomarker', inplace=True)
    output1 = output[['structure','Units', visits[0], visits[1], 'change (%)']]
    output2 = output[['structure','Units', "p-value", 'Bayes Factor', 
                      'Odds Ratio']]
    output2 = output2.astype({'Bayes Factor': 'float32'})
    output2.loc[:,'p-value'] = np.around(output2['p-value'].values, 5)
    output2.loc[:,'Bayes Factor'] = np.around(output2['Bayes Factor'].values, 2)
    output2.loc[:,'Odds Ratio'] = np.around(output2['Odds Ratio'].values, 2)
    output1.to_csv(os.path.join(src, 'stats1.csv'))
    output1.to_pickle(os.path.join(src, 'stats1.pkl'))
    output2.to_csv(os.path.join(src, 'stats2.csv'))
    output2.to_pickle(os.path.join(src, 'stats2.pkl'))



def compare_to_ref(src):

    df = pd.read_pickle(os.path.join(src, 'twoscan', 'reference.pkl'))
    df_ref = pd.read_pickle(
        os.path.join(os.getcwd(), 'tristan', 'reference.pkl'))
    
    # # This becomes obsolete afte renaming exp_med visits
    # df_ref.replace('baseline', 'control', inplace=True)
    # df_ref.replace('rifampicin', 'drug', inplace=True)

    # Add column and merge
    df_ref['source'] = 'reference'
    df['source'] = 'data'
    df = pd.concat((df_ref, df))

    # Perform t-tests and save if df output
    output = None
    for struct in df.structure.unique():
        dfs = df[df.structure==struct]
        for par in dfs.name.unique():
            dfp = dfs[dfs.name==par]
            for visit in dfp.visit.unique():
                dfv = dfp[dfp.visit==visit]
                x = dfv[dfv.source=='data'].value.values
                y = dfv[dfv.source=='reference'].value.values
                stats = pg.ttest(x, y)
                stats['structure'] = [struct]
                stats['Biomarker'] = [par]
                stats['visit'] = [visit]
                if output is None:
                    output = stats
                else:
                    output = pd.concat([output, stats])

    # Update output array
    output.reset_index(drop=True, inplace=True)
    output.drop(columns=['T', 'dof', 'alternative'], inplace=True)
    output = output[['structure','Biomarker', 'visit', 'p-val', 'CI95%',
                     'cohen-d', 'BF10', 'power']]
    output.rename(columns = {"p-val": "p-value", 'BF10': 'Bayes Factor', }, 
                  inplace=True)
    output = output.sort_values(['structure','p-value'], ascending=True)
    output.set_index('Biomarker', inplace=True)
    output = output.astype({'Bayes Factor': 'float32'})
    output.loc[:,'p-value'] = np.around(output['p-value'].values, 5)
    output.loc[:,'Bayes Factor'] = np.around(output['Bayes Factor'].values, 2)
    output.loc[:,'cohen-d'] = np.around(output['cohen-d'].values, 2)
    output.loc[:,'power'] = np.around(output['power'].values, 2)
    output = output.sort_values(['visit','structure','Biomarker'])
    output.to_csv(os.path.join(src, 'twoscan', 'stats_ref.csv'))
    output.to_pickle(os.path.join(src, 'twoscan', 'stats_ref.pkl'))



def _derive_vart_effect_sizes(output):
    # Calculate effect sizes
    for subject in output.subject.unique():
        df_subj = output[output.subject==subject]
        for structure in df_subj.structure.unique():
            df_struct = df_subj[df_subj.structure==structure]
            for tacq in df_struct.tacq.unique():
                df_tacq = df_struct[df_struct.tacq==tacq]
                for par in df_tacq.parameter.unique():
                    df_par = df_tacq[df_tacq.parameter==par]
                    if len(df_par)==2:
                        visits = df_par.visit.unique()
                        v0 = df_par[df_par.visit==visits[0]].value.values[0]
                        v1 = df_par[df_par.visit==visits[1]].value.values[0]
                        effect = 100*np.divide(v1-v0, v0)
                        name = df_par[df_par.visit==visits[0]].name.values[0]
                        unit = df_par[df_par.visit==visits[0]].unit.values[0]
                        output.loc[len(output)] = [subject, 'change (%)', structure, 
                                name, effect, unit, 0, tacq, par + ' effect size']
    return output

def derive_vart_pars(src):
    output_file = os.path.join(src, 'parameters.pkl')
    output = pd.read_pickle(output_file)
    todrop = ['BAT','BAT1','BAT2','S02b','S02','Tc']
    todrop +=['S02l','t0','t1','t2','t3','dt1','dt2',
            'khe_min','khe_max', 'khe_var',
            'kbh_min','kbh_max', 'kbh_var',
            'Khe_i','Khe_f','Kbh_i','Kbh_f', 'Th_i', 'Th_f']
    todrop += [v + ' effect size' for v in todrop]
    df = output[output.parameter.isin(todrop)]
    output.drop(index=df.index, inplace=True)
    output.reset_index(drop=True, inplace=True)
    output = _derive_vart_effect_sizes(output)
    output.to_csv(os.path.join(src, 'parameters_rep.csv'))
    output.to_pickle(os.path.join(src, 'parameters_rep.pkl'))


def first_digit(x):
    if np.isnan(x):
        return x
    return -int(math.floor(math.log10(abs(x))))

def round_sig(x, n):
    # Round to n significant digits
    if x==0:
        return x
    if np.isnan(x):
        return x
    return round(x, first_digit(x) + (n-1))

def round_to_first_digit(x):
    if np.isnan(x):
        return x
    # Round to first significant digit
    n = first_digit(x)
    y = round(x, n)
    if n<=0:
        return int(y)
    return y

def round_meas(x, xerr):
    if np.isnan(x):
        return x
    # Round measurement with known error
    if xerr==0:
        return x, xerr
    n = first_digit(xerr)
    y = round(x, n)
    yerr = round(xerr, n)
    if n<=0:
        return int(y), int(yerr)
    else:
        return y, yerr
    
def around_sig(x, n):
    return np.array([round_sig(v,n) for v in x])

def around_meas(x, xerr):
    n = len(x)
    if len(xerr) != n:
        raise ValueError(
            "The array with error values must have the same length as the "
            "array with measurements."
        )
    y = np.empty(n)
    yerr = np.empty(n)
    for i in range(n):
        y[i], yerr[i] = round_meas(x[i], xerr[i])
    return y, yerr

if __name__ == '__main__':
    assert round_to_first_digit(0.0163) == 0.02
    assert round_to_first_digit(0.163) == 0.2
    assert round_to_first_digit(1.63) == 2
    assert round_to_first_digit(16.3) == 20
    assert round_to_first_digit(163) == 200
    assert round_meas(123.654, 0.0678) == (123.65, 0.07)
    assert round_meas(123.654, 0.678) == (123.7, 0.7)
    assert round_meas(123.654, 6.78) == (124, 7)
    assert round_meas(123.654, 67.8) == (120, 70)
    assert round_meas(123.654, 678) == (100, 700)
    assert round_meas(123.654, 6780) == (0, 7000)




