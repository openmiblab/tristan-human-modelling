import math
import numpy as np

import dcmri
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig
import dcmods.pkmods as pkmods


def model(xdata,
        Te, De, ve, k_he_i, k_he_f, Th_i, Th_f, S01, S02,
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
    k_he = dcmri.interp(t, [k_he_i, k_he_f])
    # Interpolating Kbh here for consistency with original model
    Kbh = dcmri.interp(t, [1/Th_i, 1/Th_f])
    # Propagate through the gut
    ca = pkmods.flux_pfcomp(cb, Te, De, dt=dt, solver='interp')
    # Tissue concentration in the extracellular space
    Ce = ve*ca/(1-Hct)
    # Tissue concentration in the hepatocytes
    Ch = dcmri.conc_nscomp(k_he*ca, 1/Kbh, t)
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
    # Inlets
    ['Te', "Extracellular transit time", 30.0, 'sec', 0, 60, 0],
    ['De', "Extracellular dispersion", 0.85, '', 0, 1, 0],  
    # Liver tissue
    ['ve', "Liver extracellular volume fraction", 0.3, 'mL/mL', 0.01, 0.6, 0],
    ['k_he_i', "Hepatocellular uptake rate (initial)", 20/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['k_he_f', "Hepatocellular uptake rate (final)", 20/6000, 'mL/sec/mL', 0, np.inf, 0],
    ['Th_i', "Hepatocellular transit time (initial)", 30*60, 'sec', 10*60, 10*60*60, 0],
    ['Th_f', "Hepatocellular transit time (final)", 30*60, 'sec', 10*60, 10*60*60, 0],
    ['S01', "Signal amplitude S01", 1200, "a.u.", 0, np.inf, 1],
    ['S02', "Signal amplitude S02", 1200, "a.u.", 0, np.inf, 0],
]


def export_pars(pars, liver_volume, dt, t0, xdata, tR12):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters
    t = np.arange(0, max(xdata)+dt, dt)
    k_he = dcmri.interp(t, [p.k_he_i, p.k_he_f])
    Kbh = dcmri.interp(t, [1/p.Th_i, 1/p.Th_f])
    k_he_avr = np.mean(k_he)
    Kbh_avr = np.mean(Kbh)
    k_he_var = (np.amax(k_he)-np.amin(k_he))/k_he_avr
    Kbh_var = (np.amax(Kbh)-np.amin(Kbh))/Kbh_avr 
    k_bh = np.mean((1-p.ve)*Kbh)
    Th = np.mean(1/Kbh)
    df.loc['k_he'] = ["Hepatocellular uptake rate", 6000*k_he_avr, 'mL/min/100mL']
    df.loc['k_he_var'] = ["Hepatocellular uptake rate variance", 100*k_he_var, '%']
    df.loc['Kbh'] = ["Biliary tissue excretion rate", 6000*Kbh_avr, 'mL/min/100mL']
    df.loc['Kbh_var'] = ["Biliary tissue excretion rate variance", 100*Kbh_var, '%']
    df.loc['k_bh'] = ["Biliary excretion rate", 6000*k_bh, 'mL/min/100mL']
    df.loc['k_bh_i'] = ["Biliary excretion rate (initial)", 6000*(1-p.ve)/p.Th_i, 'mL/min/100mL']
    df.loc['k_bh_f'] = ["Biliary excretion rate (final)", 6000*(1-p.ve)/p.Th_f, 'mL/min/100mL']
    df.loc['Khe'] = ["Hepatocellular tissue uptake rate", 6000*k_he_avr/p.ve, 'mL/min/100mL']
    df.loc['Th'] = ["Hepatocellular mean transit time", Th/60, 'min']
    df.loc['CL_l'] = ['Liver blood clearance', 60*k_he_avr*liver_volume, 'mL/min']
    # Convert to conventional units
    df.loc['De', ['value', 'unit']] = [100*p.De, '%']
    df.loc['ve', ['value', 'unit']] = [100*p.ve, 'mL/100mL']
    df.loc['k_he_i', ['value', 'unit']] = [6000*p.k_he_i, 'mL/min/100mL']
    df.loc['k_he_f', ['value', 'unit']] = [6000*p.k_he_f, 'mL/min/100mL']
    df.loc['Th_i', ['value', 'unit']] = [p.Th_i/60, 'min']
    df.loc['Th_f', ['value', 'unit']] = [p.Th_f/60, 'min']
    # Additional exports
    # Use some notations
    t1 = xdata < tR12
    tdce1 = xdata[t1]
    t2 = xdata > tR12
    tdce2 = xdata[t2]
    df.loc['t0'] = ["Start time first acquisition", (t0+tdce1[0])/(60*60), 'hrs']
    df.loc['t1'] = ["End time first acquisition", (t0+tdce1[-1])/(60*60), 'hrs']
    df.loc['t2'] = ["Start time second acquisition", (t0+tdce2[0])/(60*60), 'hrs']
    df.loc['t3'] = ["End time second acquisition", (t0+tdce2[-1])/(60*60), 'hrs']
    df.loc['dt1'] = ["Time step first acquisition", tdce1[1]-tdce1[0], 'sec']
    df.loc['dt2'] = ["Time step second acquisition", tdce2[1]-tdce2[0], 'sec']
    return df


def plot_fit(xdata, ydata, pars, BAT1, BAT2, vars={}, xvalid=None, show=True, save=False, path=None, prefix='', xcheck=None, ycheck=None):
    clrs = ['cornflowerblue','darkblue'] 
    label = 'Liver'

    # All times
    fig.tissue_2cm(model, xdata, ydata, pars, xcheck=xcheck, ycheck=ycheck, color=clrs, label=label, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)

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