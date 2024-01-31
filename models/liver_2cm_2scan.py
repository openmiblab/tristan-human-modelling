import math
import numpy as np

import dcmri
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig
import dcmods.pkmods as pkmods


def model(xdata,
        Te, De, ve, k_he, Th, S01, S02,
        dt = 0.5,                   # Internal time resolution (sec)
        cb = None,
        Hct = 0.45,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        field_strength = 3.0,       # Field strength (T)
        R10 = 1,
        tR12 = 1,
        sample = True,
        return_conc = False,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    # Propagate through the gut
    ca = pkmods.flux_pfcomp(cb, Te, De, dt=dt, solver='interp')
    # Tissue concentration in the extracellular space
    Ce = ve*ca/(1-Hct)
    # Tissue concentration in the hepatocytes
    Ch = dcmri.conc_comp(k_he*ca, Th, t)
    if return_conc:
        return t, Ce, Ch
    # Return R
    rp = const.rp(field_strength)
    rh = const.rh(field_strength)
    R1 = R10 + rp*Ce + rh*Ch
    signal = dcmri.signal_spgress(TR, FA, R1, S01)
    t2 = (t >= tR12)
    signal[t2] = dcmri.signal_spgress(TR, FA, R1[t2], S02)
    if not sample:
        return t, signal
    return dcmri.sample(t, signal, xdata, xdata[1]-xdata[0])


def init(xdata, ydata, pars, vars, BAT1, R12):

    TR, FA = vars['TR'], vars['FA']
    R10 = vars['R10']
    tR12 = vars['tR12']

    # Use some notations
    t1 = xdata < tR12
    tdce1 = xdata[t1]
    Sdce1 = ydata[t1]
    t2 = xdata > tR12
    tdce2 = xdata[t2]
    Sdce2 = ydata[t2]

    # Estimate S01 from data
    baseline = tdce1[tdce1 <= BAT1]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, R10, 1)
    S01 = np.mean(Sdce1[:n0]) / Sref

    # Estimate S02 from data
    n0 = math.floor(60/(tdce2[1]-tdce2[0]))
    Sref = dcmri.signal_spgress(TR, FA, R12, 1)
    S02 = np.mean(Sdce2[:n0]) / Sref

    pars[-2] = S01
    pars[-1] = S02
    
    return pars, vars


PARS = [
    ['Te', "Extracellular transit time", 30.0, 'sec', 0, 60, 0],
    ['De', "Extracellular dispersion", 0.85, '', 0, 1, 0],  
    # Liver tissue
    ['ve', "Liver extracellular volume fraction", 0.3, 'mL/mL', 0.01, 0.6, 0],
    ['k_he', "Hepatocellular uptake rate", 20/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['Th', "Hepatocellular transit time", 30*60, 'sec', 10*60, 10*60*60, 0],
    ['S01', "Signal amplitude S01", 1200, "a.u.", 0, np.inf, 1],
    ['S02', "Signal amplitude S02", 1200, "a.u.", 0, np.inf, 0],
]

def export_pars(pars, liver_volume):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    df.loc['k_bh'] = ["Biliary excretion rate", 6000*(1-p.ve)/p.Th, 'mL/min/100mL']
    df.loc['Khe'] = ["Hepatocellular tissue uptake rate", 6000*p.k_he/p.ve, 'mL/min/100mL']
    df.loc['Kbh'] = ["Biliary tissue excretion rate", 6000*np.divide(1, p.Th), 'mL/min/100mL']
    df.loc['CL_l'] = ['Liver blood clearance', 60*p.k_he*liver_volume, 'mL/min']
    # Convert to conventional units
    df.loc['De', ['value', 'unit']] = [100*p.De, '%']
    df.loc['ve', ['value', 'unit']] = [100*p.ve, 'mL/100mL']
    df.loc['k_he', ['value', 'unit']] = [6000*p.k_he, 'mL/min/100mL']
    df.loc['Th', ['value', 'unit']] = [p.Th/60, 'min']
    return df


def plot_fit(xdata, ydata, pars, BAT1, BAT2, vars={}, xvalid=None, show=True, save=False, path=None, prefix='', xcheck=None, ycheck=None):
    clrs = ['cornflowerblue','darkblue'] 
    label = 'Liver'

    # All times
    fig.tissue_2cm(model, xdata, ydata, pars, color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)

    # First shot
    BAT = BAT1
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot1', xlim=[BAT-20, BAT+160], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid,  show=show, save=save, path=path, prefix=prefix)
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot1_', xlim=[BAT-20, BAT+600], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot1__', xlim=[BAT-20, BAT+1200], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)

    # Second shot
    BAT = BAT2
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot2', xlim=[BAT-20, BAT+160], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot2_', xlim=[BAT-20, BAT+600], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue_2cm(model, xdata, ydata, pars, win='shot2__', xlim=[BAT-20, BAT+1200], color=clrs, label=label, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)