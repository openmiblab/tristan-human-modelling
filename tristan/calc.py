import os
import math
import pandas as pd
import numpy as np
import pingouin as pg


def lookup(path, params, prop):
    df = pd.read_csv(os.path.join(path, 'output_data_dict.csv'))
    df.set_index('parameter', inplace=True)
    vals = []
    for p in params:
        if ' effect size' in p:
            if prop=='unit':
                vals.append('%')
            else:
                vals.append(df.at[p[:-12], prop])
        else:
            vals.append(df.at[p, prop])
    return vals


def _derive_effect_sizes(output):
    effect_size = []
    for subject in output.subject.unique():
        df_subj = output[output.subject==subject]
        for par in df_subj.parameter.unique():
            df_par = df_subj[df_subj.parameter==par]
            if len(df_par)==2:
                visits = df_par.visit.unique()
                v0 = df_par[df_par.visit==visits[0]].value.values[0]
                v1 = df_par[df_par.visit==visits[1]].value.values[0]
                effect = 100*np.divide(v1-v0, v0)
                effect = [subject, par, effect]
                effect_size.append(effect)
    return pd.DataFrame(
        data=effect_size, 
        columns=['subject','parameter','value'],
    )

def _derive_vart_effect_sizes(output):
    effect_size = []
    for subject in output.subject.unique():
        df_subj = output[output.subject==subject]
        for tacq in df_subj.tacq.unique():
            df_tacq = df_subj[df_subj.tacq==tacq]
            for par in df_tacq.parameter.unique():
                df_par = df_tacq[df_tacq.parameter==par]
                if len(df_par)==2:
                    visits = df_par.visit.unique()
                    v0 = df_par[df_par.visit==visits[0]].value.values[0]
                    v1 = df_par[df_par.visit==visits[1]].value.values[0]
                    effect = 100*np.divide(v1-v0, v0)
                    effect = [subject, par, effect, tacq]
                    effect_size.append(effect)
    return pd.DataFrame(
        data=effect_size, 
        columns=['subject','parameter','value', 'tacq'],
    )


def derive_pars(src, ref=False):
    output_file = os.path.join(src, 'output_data.csv')
    output = pd.read_csv(output_file)
    todrop = [
        'BAT','BAT1','BAT2','S02a','S02',
        'S02l','t0','t1','t2','t3','dt1','dt2',
        'khe_min','khe_max', 'khe_var',
        'kbh_min','kbh_max', 'kbh_var',
        'Khe_i','Khe_f','Kbh_i','Kbh_f', 'Th_i', 'Th_f',
    ]
    df = output[output.parameter.isin(todrop)]
    output.drop(index=df.index, inplace=True)
    output.reset_index(drop=True, inplace=True)
    path = os.path.join(src, 'Analysis')
    if not os.path.exists(path):
        os.makedirs(path)
    if ref:
        output.to_csv(os.path.join(path, 'reference.csv'), index=False)
    else:
        output.to_csv(os.path.join(path, 'parameters_rep.csv'), index=False)
        effect_size = _derive_effect_sizes(output)
        effect_size.to_csv(os.path.join(path, 'effect_size.csv'), index=False)


def derive_vart_pars(src):
    output_file = os.path.join(src, 'output_data.csv')
    output = pd.read_csv(output_file)
    todrop =[
        'BAT','BAT1','BAT2','S02b','S02','Tc'
        'S02l',
        't0','t1','t2','t3','dt1','dt2',
        'khe_min','khe_max', 'khe_var',
        'kbh_min','kbh_max', 'kbh_var',
        'Khe_i','Khe_f','Kbh_i','Kbh_f', 'Th_i', 'Th_f',
    ]
    df = output[output.parameter.isin(todrop)]
    output.drop(index=df.index, inplace=True)
    output.reset_index(drop=True, inplace=True)
    path = os.path.join(src, 'Analysis')
    if not os.path.exists(path):
        os.makedirs(path)
    output.to_csv(os.path.join(path, 'parameters_rep.csv'), index=False)
    effect_size = _derive_vart_effect_sizes(output)
    effect_size.to_csv(os.path.join(path, 'effect_size.csv'), index=False)
    



def desc_stats(src):
    path = os.path.join(src, 'Analysis')
    if not os.path.exists(path):
        os.makedirs(path)
    output_file = os.path.join(path, 'parameters_rep.csv')
    output = pd.read_csv(output_file)
    visits = output.visit.unique()
    v0 = output[output.visit==visits[0]]
    v0 = pd.pivot_table(v0, values='value', index='subject', 
                        columns='parameter')
    v1 = output[output.visit==visits[1]]
    v1 = pd.pivot_table(v1, values='value', index='subject', 
                        columns='parameter')
    
    effect_file = os.path.join(path, 'effect_size.csv')
    ef = pd.read_csv(effect_file)
    ef = pd.pivot_table(ef, values='value', index='subject', 
                        columns='parameter')
    v0.to_csv(os.path.join(path, '_table_'+visits[0]+'.csv'))
    v1.to_csv(os.path.join(path, '_table_'+visits[1]+'.csv'))
    ef.to_csv(os.path.join(path, '_table_effect.csv'))
    # Calculate stats
    bstats = v0.describe()
    bstats = bstats[['khe', 'kbh']].round(2)
    bstats = bstats.rename(columns={
        "khe": "khe "+visits[0]+" (mL/min/100cm3)", 
        "kbh": "kbh "+visits[0]+" (mL/min/100cm3)"})
    rstats = v1.describe()
    rstats = rstats[['khe', 'kbh']].round(2)
    rstats = rstats.rename(columns={
        "khe": "khe "+visits[1]+" (mL/min/100cm3)", 
        "kbh": "kbh "+visits[1]+" (mL/min/100cm3)"})
    estats = ef.describe()
    estats = estats[['khe', 'kbh']].round(1)
    estats = estats.rename(columns={
        "khe": "khe effect size (%)", 
        "kbh": "kbh effect size (%)"})
    stats = pd.concat([estats.T, bstats.T, rstats.T])
    stats=stats.reindex([
        "khe effect size (%)",
        "khe "+visits[0]+" (mL/min/100cm3)",
        "khe "+visits[1]+" (mL/min/100cm3)",
        "kbh effect size (%)",
        "kbh "+visits[0]+" (mL/min/100cm3)",
        "kbh "+visits[1]+" (mL/min/100cm3)",
        ])
    stats.to_csv(os.path.join(path, '_table_k_stats.csv'))



def ttest(src):
    
    path = os.path.join(src, 'Analysis')
    if not os.path.exists(path):
        os.makedirs(path)

    df = pd.read_csv(os.path.join(path, 'parameters_rep.csv'))
    df['group'] = lookup(src, df.parameter.values, 'group')
    #df['name'] = lookup(src, df.parameter.values, 'name')
    df['unit'] = lookup(src, df.parameter.values, 'unit') 
    visits = df.visit.unique()
    eff = pd.read_csv(os.path.join(path, 'effect_size.csv'))
    eff['visit'] = 'change (%)'
    eff['group'] = lookup(src, eff.parameter.values, 'group')
    #eff['name'] = lookup(src, eff.parameter.values, 'name')
    eff['unit'] = lookup(src, eff.parameter.values, 'unit') 
    
    df = df[['subject','visit','group','parameter','value','unit']]
    eff = eff[['subject','visit','group','parameter','value','unit']]
    df = pd.concat([df, eff])

    # Get mean, sdev and count of all variables
    avr = pd.pivot_table(df, values='value', columns='visit', 
                         index=['group','parameter','unit'], aggfunc='mean')
    std = pd.pivot_table(df, values='value', columns='visit', 
                         index=['group','parameter','unit'], aggfunc='std')
    cnt = pd.pivot_table(df, values='value', columns='visit', 
                         index=['group','parameter','unit'], aggfunc='count')

    avr = avr.sort_values(['group','parameter']) 
    std = std.sort_values(['group','parameter'])
    cnt = cnt.sort_values(['group','parameter'])

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
    for struct in df.group.unique():
        dfs = df[df.group==struct]
        for par in dfs.parameter.unique():
            dfp = dfs[dfs.parameter==par]
            stats = pg.pairwise_tests(data=dfp, 
                    dv='value', within='visit', subject='subject', 
                    return_desc=False, effsize='odds-ratio')
            stats['group'] = [struct]
            stats['parameter'] = [par]
            if output is None:
                output = stats
            else:
                output = pd.concat([output, stats])
    output = output.sort_values(['group','parameter'])

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
    output['name'] = lookup(src, output.parameter.values, 'name')
    output = output[['group','name', 'unit', visits[0], visits[1], 
                     'change (%)', 'p-unc','BF10', 'odds-ratio']]
    output.rename(columns={"name": "Biomarker", "unit": "Units", 
                           "p-unc": "p-value", 'BF10': 'Bayes Factor', 
                           'odds-ratio': 'Odds Ratio'}, inplace=True)
    output = output.sort_values(['group','p-value'], ascending=True)
    output.set_index('Biomarker', inplace=True)
    output1 = output[['group','Units', visits[0], visits[1], 'change (%)']]
    output2 = output[['group','Units', "p-value", 'Bayes Factor', 
                      'Odds Ratio']]
    output2 = output2.astype({'Bayes Factor': 'float32'})
    output2.loc[:,'p-value'] = np.around(output2['p-value'].values, 5)
    output2.loc[:,'Bayes Factor'] = np.around(output2['Bayes Factor'].values, 2)
    output2.loc[:,'Odds Ratio'] = np.around(output2['Odds Ratio'].values, 2)
    output1.to_csv(os.path.join(path, 'stats1.csv'))
    output2.to_csv(os.path.join(path, 'stats2.csv'))



def compare_to_ref(src):

    df = pd.read_csv(os.path.join(src, 'output_data.csv'))
    df_ref = pd.read_csv(
        os.path.join(os.getcwd(), 'tristan', 'reference.csv')
    )

    # Add column and merge
    df['source'] = 'data'
    df_ref['source'] = 'reference'
    df = pd.concat((df_ref, df))
    df['group'] = lookup(src, df.parameter.values, 'group')

    # Perform t-tests and save if df output
    output = None
    for struct in df.group.unique():
        dfs = df[df.group==struct]
        for par in dfs.parameter.unique():
            dfp = dfs[dfs.parameter==par]
            for visit in dfp.visit.unique():
                dfv = dfp[dfp.visit==visit]
                if dfv.source.unique().size==1:
                    continue
                stats = pg.pairwise_tests(data=dfv, 
                        dv='value', within='source', subject='subject', 
                        return_desc=False, effsize='odds-ratio')
                stats['group'] = [struct]
                stats['parameter'] = [par]
                stats['visit'] = [visit]
                if output is None:
                    output = stats
                else:
                    output = pd.concat([output, stats])
    output = output.sort_values(['visit','group','parameter'])

    # Update output array
    output.reset_index(drop=True, inplace=True)
    output.drop(columns=['Contrast', 'A', 'B', 'Paired', 'Parametric', 'T', 
                         'dof', 'alternative'], inplace=True)
    output['name'] = lookup(src, output.parameter.values, 'name')
    output['unit'] = lookup(src, output.parameter.values, 'unit') 
    output = output[['group','name', 'unit', 'visit', 
                     'p-unc', 'BF10', 'odds-ratio']]
    output.rename(columns={"name": "Biomarker", "unit": "Units", 
                           "p-unc": "p-value", 'BF10': 'Bayes Factor', 
                           'odds-ratio': 'Odds Ratio'}, inplace=True)
    output = output.sort_values(['group','p-value'], ascending=True)
    output.set_index('Biomarker', inplace=True)
    output = output.astype({'Bayes Factor': 'float32'})
    output.loc[:,'p-value'] = np.around(output['p-value'].values, 5)
    output.loc[:,'Bayes Factor'] = np.around(output['Bayes Factor'].values, 2)
    output.loc[:,'Odds Ratio'] = np.around(output['Odds Ratio'].values, 2)
    output.to_csv(os.path.join(src, 'Analysis', 'stats_ref.csv'))






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




