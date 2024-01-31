import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        vb, r, FF, Tp,
        R10k = 1,
        S0 = 1,                     # Baseline signal (a.u.)
        dt = 0.5,                   # Internal time resolution (sec)
        cb = None,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        field_strength = 3.0,       # Field strength (T)
        agent = 'Dotarem',
        Hct = 0.45, 
        return_conc = False,
        sample = True,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    vp = (1-Hct)*vb
    Fp = vp/Tp/(1+r*FF)
    Ft = FF*Fp
    vt = 1-vb
    Fu = (1-r)*Ft
    Tt = vt/Fu
    ca = cb/(1-Hct)
    Cp = dcmri.conc_comp(Fp*ca, Tp, t)
    cp = Cp/vp
    Ct = dcmri.conc_comp(Ft*cp, Tt, t)
    if return_conc:
        return t, Cp, Ct
    # Return R
    rp = const.rp(field_strength, agent)
    R1k = R10k + rp*Cp + rp*Ct
    signal = dcmri.signal_spgress(TR, FA, R1k, S0)
    if not sample:
        return t, signal
    return dcmri.sample(t, signal, xdata, xdata[1]-xdata[0])


def init(xdata, ydata, vars, BAT):
    # Estimate S0 from data
    baseline = xdata[xdata <= BAT]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(vars['TR'], vars['FA'], vars['R10k'], 1)
    S0 = np.mean(ydata[:n0]) / Sref
    vars['S0'] = S0
    return vars


PARS = [
    ['vb', "Blood volume", 0.3, 'mL/mL', 0.01, 0.6, 0],
    ['Tp', "Plasma transit time", 5, 'sec', 0, 8, 0],
    ['FF', "Filtration Fraction", 0.1, '', 0.01, 0.5, 0],
    ['r', "Reabsorption fraction", 0.95, '', 0.5, 1.0, 0],
]


def export_pars(pars, vars):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    vp = (1-vars['Hct'])*p.vb
    Fp = vp/p.Tp/(1+p.r*p.FF)
    Ft = p.FF*Fp
    vt = 1-p.vb
    Fu = (1-p.r)*Ft
    Fb = Fp/(1-vars['Hct'])
    df.loc['Fb'] = ["Renal blood flow", 6000*Fb, 'mL/min/100mL']
    df.loc['Ft'] = ["Tubular flow", 6000*Ft, 'mL/min/100mL']
    df.loc['Fu'] = ["Urine flow", 6000*Fu, 'mL/min/100mL']    
    df.loc['vt'] = ["Tubular volume fraction", 100*vt, 'mL/100mL']
    df.loc['Tt'] = ["Tubular transit time", vt/Fu/60, 'min']
    # Convert units
    df.loc['vb', ['value', 'unit']] = [100*p.vb, 'mL/100mL']
    df.loc['FF', ['value', 'unit']] = [100*p.FF, '%']
    df.loc['r', ['value', 'unit']] = [100*p.r, '%']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
