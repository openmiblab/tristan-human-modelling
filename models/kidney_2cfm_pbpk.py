import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        Ta, FF_k, F_k, Tp, Tt,
        R10k = 1,
        S0 = 1,                     # Baseline signal (a.u.)
        dt = 0.5,                   # Internal time resolution (sec)
        J_aorta = None,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        field_strength = 3.0,       # Field strength (T)
        agent = 'Dotarem',
        kidney_volume = None, 
        return_conc = False,
        sample = True,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    E_k = F_k/(1+F_k)
    J_kidneys = FF_k*J_aorta
    J_kidneys = dcmri.flux_plug(J_kidneys, Ta, t)
    Np = dcmri.conc_comp(J_kidneys, Tp, t)
    Nt = dcmri.conc_comp(E_k*Np/Tp, Tt, t)
    Cp = 1000*Np/kidney_volume # mM
    Ct = 1000*Nt/kidney_volume # mM
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
    ['Ta',"Arterial delay time", 1.0, "sec", 0.0, 3.0, 0],
    ['FF_k', "Kidney flow fraction", 0.2, "", 0.01, 1, 0],
    ['F_k', "Kidney filtration fraction", 0.10, "", 0.0, 1.0, 0],
    ['Tp', "Vascular mean transit time", 3.0, "sec", 1.0, 10.0, 0],
    ['Tt', "Tubular mean transit time", 90.0, "sec", 1.0, np.inf, 0],
]


def export_pars(pars, kidney_volume, CO, Hct):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    E_k = np.divide(p.F_k, 1+p.F_k)
    df.loc['E_k'] = ["Kidney extraction fraction", 100*E_k, '%']
    df.loc['GFR'] = ["Glomerular Filtration Rate", 60*p.F_k*p.FF_k*CO*(1-Hct), 'mL/min']  
    df.loc['RBF'] = ["Renal blood flow", 60*p.FF_k*CO, 'mL/min']
    df.loc['RP'] = ["Renal perfusion", 6000*p.FF_k*CO/kidney_volume, 'mL/min/100mL']
    df.loc['RBV'] = ["Renal blood volume", 100*p.FF_k*CO*p.Tp/kidney_volume, 'mL/100mL']
    # Convert to conventional units
    df.loc['FF_k', ['value', 'unit']] = [100*p.FF_k, '%']
    df.loc['F_k', ['value', 'unit']] = [100*p.F_k, '%']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
