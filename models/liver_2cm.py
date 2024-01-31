import numpy as np

import dcmri
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig
import dcmods.pkmods as pkmods


def model(xdata,
        Te, De, ve, k_he, Th, S0,
        dt = 0.5,                   # Internal time resolution (sec)
        cb = None,
        Hct = 0.45,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        R10 = 1,
        field_strength = 3.0,       # Field strength (T)
        sample = True,
        return_conc = False,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    # Propagate through the gut
    ca = pkmods.flux_pfcomp(cb, Te, De, dt=dt, solver='interp')
    #ca = dcmri.flux_comp(cb, Te, t)
    #ca = dcmri.flux_plug(ca, Td, t, solver='interp')
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
    signal = dcmri.signal_spgress(TR, FA, R1, S0)
    if not sample:
        return t, signal
    return dcmri.sample(t, signal, xdata, xdata[1]-xdata[0])


def init(xdata, ydata, pars, vars, BAT):
    R10 = vars['R10']
    tdce, Sdce = xdata, ydata
    TR, FA = vars['TR'], vars['FA']
    
    # Estimate S0 from data
    baseline = tdce[tdce <= BAT]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, R10, 1)
    S0 = np.mean(Sdce[:n0]) / Sref
    pars[-1] = S0
    
    return pars, vars


PARS = [
    # Inlets
    ['Te', "Extracellular transit time", 30.0, 'sec', 0.1, 60, 0], 
    ['De', "Extracellular dispersion", 0.85, '', 0, 1, 0],
    # Liver tissue
    ['ve', "Liver extracellular volume fraction", 0.3, 'mL/mL', 0.01, 0.6, 0],
    ['k_he', "Hepatocellular uptake rate", 20/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['Th', "Hepatocellular transit time", 30*60, 'sec', 10*60, 10*60*60, 0],
    ['S0', "Signal amplitude S0", 1200, "a.u.", 0, np.inf, 1],
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


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix='', xcheck=None, ycheck=None):
    clrs = ['cornflowerblue','darkblue'] 
    label = 'Liver'
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, xcheck=xcheck, ycheck=ycheck, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, xcheck=xcheck, ycheck=ycheck, win='win', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, xcheck=xcheck, ycheck=ycheck, win='win_', xlim=[BAT-20, BAT+600], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, xcheck=xcheck, ycheck=ycheck, win='win__', xlim=[BAT-20, BAT+1200], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)