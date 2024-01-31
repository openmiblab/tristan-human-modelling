import numpy as np

import dcmri

import dcmods.const as const
import dcmods.pkmods as pkmods
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        BAT1, BAT2, CO, T_hl, D_hl, E_o, Tp_o, Te_o, E_b, S01, S02,
        dt = 0.5,                   # Internal time resolution (sec)
        weight = 70.0,              # Patient weight in kg
        conc = 0.25,                # mmol/mL (https://www.bayer.com/sites/default/files/2020-11/primovist-pm-en.pdf)
        dose = [0.025, 0.025],      # mL per kg bodyweight (quarter dose)
        rate = 1,                   # Injection rate (mL/sec)
        dose_tolerance = 0.1,
        field_strength = 3.0,       # Field strength (T)
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        R10 = 1.0,                  # Precontrast relaxation rate (1/sec)
        tR12 = 1,                  # 
        sample = True,              # If false return the pseudocontinuous signal
        return_conc = False,
    ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    Ji = dcmri.inj_step(t, 
        weight, conc, dose[0], rate, BAT1,
        dose[1], BAT2)
    Jb = pkmods.flux_wb_pcc(Ji,
            T_hl, D_hl,
            E_o, Tp_o, Te_o,
            E_b, 
            dt=dt, tol=dose_tolerance)
    cb = Jb*1000/CO  # (mM)
    if return_conc:
        return t, cb
    rp = const.rp(field_strength)
    R1 = R10 + rp*cb
    signal = dcmri.signal_spgress(TR, FA, R1, S01)
    t2 = (t >= tR12)
    signal[t2] = dcmri.signal_spgress(TR, FA, R1[t2], S02)
    if not sample:
        return t, signal
    return dcmri.sample(t, signal, xdata, xdata[1]-xdata[0])


def init(xdata, ydata, pars, vars, R12):

    TR, FA = vars['TR'], vars['FA']
    R10 = vars['R10']
    tR12 = vars['tR12']

    # Use some notations
    k = xdata < tR12
    tdce1 = xdata[k]
    Sdce1 = ydata[k]
    k = xdata > tR12
    tdce2 = xdata[k]
    Sdce2 = ydata[k]
    T, D = pars[3], pars[4]

    # Estimate BAT1 and S01 from data
    BAT1 = tdce1[np.argmax(Sdce1)] - (1-D)*T
    baseline = tdce1[tdce1 <= BAT1-20]
    n0 = baseline.size
    if n0==0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, R10, 1)
    S01 = np.mean(Sdce1[:n0]) / Sref

    # Estimate BAT2 and S02 from data
    BAT2 = tdce2[np.argmax(Sdce2)] - (1-D)*T
    baseline = tdce2[tdce2 <= BAT2-20]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, R12, 1)
    S02 = np.mean(Sdce2[:n0]) / Sref

    # Format data in standard form
    pars[0] = BAT1
    pars[1] = BAT2
    pars[-2] = S01
    pars[-1] = S02

    # Adjust dt if needed
    dose = min(vars['dose'])
    dose = dose if dose>0 else max(vars['dose'])
    duration = vars['weight']*dose/vars['rate']
    if duration != 0:
        if vars['dt'] > duration/5:
            vars['dt'] = duration/5
    
    return pars, vars


PARS = [
    ['BAT1', "Bolus arrival time 1", 60, "sec", 0, np.inf, 0],
    ['BAT2', "Bolus arrival time 2", 1200, "sec", 0, np.inf, 0],
    ['CO', "Cardiac output", 100.0, "mL/sec", 0, np.inf, 0], # 6 L/min = 100 mL/sec
    ['T_hl', "Heart-lung mean transit time", 10.0, "sec", 0, 30, 0],
    ['D_hl', "Heart-lung transit time dispersion", 0.2, "", 0.05, 0.95, 0],
    ['E_o', "Extracellular extraction fraction", 0.15, "", 0, 0.50, 0],
    ['Tp_o', "Organs mean transit time", 20.0, "sec", 0, 60, 0],
    ['Te_o', "Extracellular mean transit time", 120.0, "sec", 0, 800.0, 0],
    ['E_b',"Body extraction fraction", 0.05,"", 0.01, 0.15, 0],
    ['S01', "Signal amplitude S01", 1200, "a.u.", 0, np.inf, 1],
    ['S02', "Signal amplitude S02", 1200, "a.u.", 0, np.inf, 0],
]

def export_pars(pars):
    df = tools.df_export(PARS, pars)
    v = df.value 
    # Add derived parameters 
    df.loc['Tc'] = ["Mean circulation time", v.Tp_o+v.T_hl, 'sec'] 
    # Export units
    df.loc['E_o', ['value', 'unit']] = [100*v.E_o, '%']
    df.loc['E_b', ['value', 'unit']] = [100*v.E_b, '%']
    return df


def plot_fit(xdata, ydata, pars, vars={}, xvalid=None, show=True, save=False, path=None, prefix='', xcheck=None, ycheck=None):
    clr = ['lightcoral','darkred']
    lbl = 'Aorta'
    fig.tissue(model, xdata, ydata, pars, color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    BAT = pars[0]
    fig.tissue(model, xdata, ydata, pars, win='shot1', xlim=[BAT-20, BAT+160], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue(model, xdata, ydata, pars, win='shot1_', xlim=[BAT-20, BAT+600], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue(model, xdata, ydata, pars, win='shot1__', xlim=[BAT-20, BAT+1200], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    BAT = pars[1]
    fig.tissue(model, xdata, ydata, pars, win='shot2', xlim=[BAT-20, BAT+160], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue(model, xdata, ydata, pars, win='shot2_', xlim=[BAT-20, BAT+600], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    fig.tissue(model, xdata, ydata, pars, win='shot2__', xlim=[BAT-20, BAT+1200], color=clr, label=lbl, xcheck=xcheck, ycheck=ycheck, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
