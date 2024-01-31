import numpy as np

import dcmri

import dcmods.const as const
import dcmods.pkmods as pkmods
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        BAT, CO, T_hl, D_hl, E_o, Tp_o, Te_o, E_b,
        S0 = 1,                     # Baseline signal (a.u.)
        dt = 0.5,                   # Internal time resolution (sec)
        weight = 70.0,              # Patient weight in kg
        conc = 0.25,                # mmol/mL (https://www.bayer.com/sites/default/files/2020-11/primovist-pm-en.pdf)
        dose = 0.025,               # mL per kg bodyweight (quarter dose)
        rate = 1,                   # Injection rate (mL/sec)
        dose_tolerance = 0.1,
        field_strength = 3.0,       # Field strength (T)
        agent = 'Dotarem',
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        R10 = 1.0,                  # Precontrast relaxation rate (1/sec)
        sample = True,              # If false return the pseudocontinuous signal
        return_conc = False,
    ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    Ji = dcmri.inj_step(t, 
        weight, conc, dose, rate, BAT)
    Jb = pkmods.flux_wb_chc(Ji,
            T_hl, D_hl,
            E_o, Tp_o, Te_o,
            E_b, 
            dt=dt, tol=dose_tolerance)
    cb = Jb*1000/CO  # (mM)
    if return_conc:
        return t, cb
    rp = const.rp(field_strength, agent)
    R1 = R10 + rp*cb
    signal = dcmri.signal_spgress(TR, FA, R1, S0)
    if not sample:
        return t, signal
    return dcmri.sample(t, signal, xdata, xdata[1]-xdata[0])


PARS = [
    ['BAT', "Bolus arrival time", 60, "sec", 0, np.inf, 0],
    ['CO', "Cardiac output", 100.0, "mL/sec", 0, np.inf, 0], # 6 L/min = 100 mL/sec
    ['T_hl', "Heart-lung mean transit time", 10.0, "sec", 0, 30, 0],
    ['D_hl', "Heart-lung transit time dispersion", 0.2, "", 0.05, 0.95, 0],
    ['E_o', "Extracellular extraction fraction", 0.15, "", 0, 0.50, 0],
    ['Tp_o', "Organs mean transit time", 20.0, "sec", 0, 60, 0],
    ['Te_o', "Extracellular mean transit time", 120.0, "sec", 0, 800.0, 0],
    ['E_b',"Body extraction fraction", 0.05,"", 0.01, 0.15, 0],
]


def init(xdata, ydata, pars, vars):

    TR, FA = vars['TR'], vars['FA']

    # Estimate BAT and S0 from data
    BAT = xdata[np.argmax(ydata)]
    baseline = xdata[xdata <= BAT-20]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, vars['R10'], 1)
    S0 = np.mean(ydata[:n0]) / Sref

    # Adjust dt if needed
    duration = vars['weight']*vars['dose']/vars['rate']
    if duration != 0:
        if vars['dt'] > duration/5:
            vars['dt'] = duration/5

    # Initialize variables
    pars[0] = BAT
    vars['S0'] = S0 
    
    return ydata, pars, vars


def export_pars(pars):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    df.loc['Tc'] = ["Mean circulation time", p.Tp_o+p.T_hl, 'sec'] 
    # Convert to conventional units
    df.loc['E_o', ['value', 'unit']] = [100*p.E_o, '%']
    df.loc['E_b', ['value', 'unit']] = [100*p.E_b, '%']
    return df


def plot_fit(xdata, ydata, pars, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clr = ['lightcoral','darkred']
    lbl = 'Aorta'
    fig.tissue(model, xdata, ydata, pars, color=clr, label=lbl, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)
    BAT = pars[0]
    fig.tissue(model, xdata, ydata, pars, win='pass1', xlim=[BAT-20, BAT+160], color=clr, label=lbl, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix)