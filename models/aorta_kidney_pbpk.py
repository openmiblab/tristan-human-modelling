import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        BAT, CO, T_hl, D_hl, E_o, Tp_o, Te_o, E_l,
        Ta, FF_k, F_k, Tp, dr, ut0, ut1, ut2, r0, r1, r2,
        R10b = 1,
        R10k = 1,
        S0k = 1,                     # Baseline signal (a.u.)
        S0b = 1,
        dt = 0.5,                   # Internal time resolution (sec)
        weight = 70.0,              # Patient weight in kg
        conc = 0.25,                # mmol/mL (https://www.bayer.com/sites/default/files/2020-11/primovist-pm-en.pdf)
        dose = 0.025,               # mL per kg bodyweight (quarter dose)
        rate = 1,                   # Injection rate (mL/sec)
        dose_tolerance = 0.1,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        field_strength = 3.0,       # Field strength (T)
        agent = 'Dotarem',
        kidney_volume = None, 
        return_conc = False,
        sample = True,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    # Aorta
    Ji = dcmri.inj_step(t, # mmol/sec
        weight, conc, dose, rate, BAT)
    J_aorta = pkmods.prop_body3_kidneys(t, Ji,
        T_hl, D_hl,
        E_o, Tp_o, Te_o,
        E_l,
        FF_k, E_k, 1/Tp, Ta,
        tol=dose_tolerance)
    cb = 1000*J_aorta/CO # mM
    # Kidney
    E_k = F_k/(1+F_k)
    vel = [ut0, ut1, ut2]
    dif = [dr]
    reabs = [r0, r1, r2]
    J_kidneys = FF_k*J_aorta
    J_kidneys = dcmri.flux_plug(J_kidneys, Ta, t)
    Np = dcmri.conc_nscomp(J_kidneys, Tp, t)
    Nt = pkmods.conc_neph(t, E_k*Np/Tp, vel, dif, reabs, n=50)
    Cp = 1000*Np/kidney_volume # mM
    Ct = 1000*Nt/kidney_volume # mM
    if return_conc:
        return t, cb, Cp, Ct
    # Return R
    rp = const.rp(field_strength, agent)
    R1k = R10k + rp*Cp + rp*Ct
    R1b = R10b + rp*cb
    Sb = dcmri.signal_spgress(TR, FA, R1b, S0b)
    Sk = dcmri.signal_spgress(TR, FA, R1k, S0k)
    if not sample:
        return t, Sb, Sk
    Sb = dcmri.sample(t, Sb, xdata, xdata[1]-xdata[0])
    Sk = dcmri.sample(t, Sk, xdata, xdata[1]-xdata[0])
    return np.concatenate([Sb, Sk])


def init(xdata, ydata, pars, vars):

    # Adjust dt if needed
    duration = vars['weight']*vars['dose']/vars['rate']
    if duration != 0:
        if vars['dt'] > duration/5:
            vars['dt'] = duration/5

    n = xdata.size/2
    x = xdata[:n]
    yb = ydata[:n]
    yk = ydata[n:]

    TR, FA = vars['TR'], vars['FA']

    # Estimate BAT and aorta S0 from data
    BAT = x[np.argmax(yb)]
    baseline = x[x <= BAT-20]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, vars['R10b'], 1)
    S0b = np.mean(yb[:n0]) / Sref

    # Estimate kidney S0 from data
    Sref = dcmri.signal_spgress(TR, FA, vars['R10k'], 1)
    S0k = np.mean(yk[:n0]) / Sref

    pars[0] = BAT
    vars['S0b'] = S0b
    vars['S0k'] = S0k
    return pars, vars


PARS = [
    # Aorta
    ['BAT', "Bolus arrival time", 60, "sec", 0, np.inf, 0],
    ['CO', "Cardiac output", 100.0, "mL/sec", 0, np.inf, 0], # 6 L/min = 100 mL/sec
    ['T_hl', "Heart-lung mean transit time", 10.0, "sec", 0, 30, 0],
    ['D_hl', "Heart-lung transit time dispersion", 0.2, "", 0.05, 0.95, 0],
    ['E_o', "Extracellular extraction fraction", 0.15, "", 0, 0.50, 0],
    ['Tp_o', "Organs mean transit time", 20.0, "sec", 0, 60, 0],
    ['Te_o', "Extracellular mean transit time", 120.0, "sec", 0, 800.0, 0],
    ['E_l',"Liver extraction fraction", 0.05,"", 0.01, 0.15, 0],
    # Kidney
    ['Ta',"Arterial delay time", 1.0,"sec", 0.0, 3.0, 0],
    ['FF_k', "Kidney flow fraction", 0.2, "", 0.01, 1, 0],
    ['F_k', "Kidney filtration fraction", 0.10, "", 0.0, 1.0, 0],
    ['Tp', "Plasma mean transit time", 3.0,"sec", 1.0, 10.0, 0],
    ['dr', "Diffusion rate", 1e-4, '1/sec', 0, 1e-2, 0],
    ['ut0', "Tubular velocity", 1/100., '1/sec', 0, 1/5, 0],
    ['ut1', "Tubular velocity", 1/100., '1/sec', 0, 1/5, 0],
    ['ut2', "Tubular velocity", 1/100., '1/sec', 0, 1/5, 0],
    ['r0', "Reabsorption rate", -np.log(1-0.95), '', 0, np.inf, 0],
    ['r1', "Reabsorption rate", -np.log(1-0.95), '', 0, np.inf, 0],
    ['r2', "Reabsorption rate", -np.log(1-0.95), '', 0, np.inf, 0],
]


def export_pars(pars, kidney_volume, tmax, CO, Hct):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    E_k = np.divide(p.F_k, 1+p.F_k)
    df.loc['E_k'] = ["Kidney extraction fraction", 100*E_k, '%']
    df.loc['GFR'] = ["Glomerular Filtration Rate", 60*p.F_k*p.FF_k*CO*(1-Hct), 'mL/min']  
    df.loc['RBF'] = ["Renal blood flow", 60*p.FF_k*CO, 'mL/min']
    df.loc['RP'] = ["Renal perfusion", 6000*p.FF_k*CO/kidney_volume, 'mL/min/100mL']
    df.loc['RBV'] = ["Renal blood volume", 100*p.FF_k*CO*p.Tv/kidney_volume, 'mL/100mL']
    # Convert to conventional units
    df.loc['E_o', ['value', 'unit']] = [100*p.E_o, '%']
    df.loc['E_l', ['value', 'unit']] = [100*p.E_l, '%']
    df.loc['FF_k', ['value', 'unit']] = [100*p.FF_k, '%']
    df.loc['F_k', ['value', 'unit']] = [100*p.F_k, '%']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    # This needs a new function plotting aorta and kidney
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
