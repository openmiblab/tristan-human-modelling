import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        Tp, Fp, Ft, Tt,
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
    vp = Tp*(Fp+Ft)
    ca = cb/(1-Hct)
    cp = dcmri.conc_comp(Fp*ca, Tp, t)
    Cp = vp*cp
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
    ['Fp', "Plasma flow", 200/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['Tp', "Plasma mean transit time", 5, 'sec', 0.0, 8, 0],
    ['Ft', "Tubular flow", 30/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['Tt', "Tubular mean transit time", 120, 'sec', 1, np.inf, 0],
]


def export_pars(pars, vars):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    df.loc['vp'] = ["Plasma volume fraction", 100*p.Tp*(p.Fp+p.Ft), 'mL/100mL']         
    df.loc['FF'] = ["Filtration Fraction", 100*p.Ft/p.Fp, '%']
    # Convert to conventional units
    df.loc['Fp', ['value', 'unit']] = [6000*p.Fp, 'mL/min/100mL']
    df.loc['Ft', ['value', 'unit']] = [6000*p.Ft, 'mL/min/100mL']
    df.loc['Tt', ['value', 'unit']] = [p.Tt/60, 'min']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
