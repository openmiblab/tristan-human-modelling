import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig

#liver: 5
#aorta: 8
#combined: 14 (+1)
def model(xdata,
        BAT, CO, T_hl, D_hl, E_o, Tp_o, Te_o, 
        FF_l, E_l, E_k, 
        Te, De, Tel, Th,
        R10b = 1,
        R10l = 1,
        S0l = 1,                    # Baseline signal (a.u.)
        S0b = 1,
        dt = 0.5,                   # Internal time resolution (sec)
        weight = 70.0,              # Patient weight in kg
        conc = 0.25,                # mmol/mL (https://www.bayer.com/sites/default/files/2020-11/primovist-pm-en.pdf)
        dose = 0.025,               # mL per kg bodyweight (quarter dose)
        rate = 1,                   # Injection rate (mL/sec)
        dose_tolerance = 0.1,
        TR = 3.71/1000.0,           # Repetition time (sec)
        FA = 15.0,                  # Nominal flip angle (degrees)
        Hct = 0.45, 
        field_strength = 3.0,       # Field strength (T)
        liver_volume = None, 
        return_conc = False,
        sample = True,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    # Aorta
    Ji = dcmri.inj_step(t, weight, conc, dose, rate, BAT) # mmol/sec
    J_aorta = pkmods.flux_wb_liver(Ji,
        T_hl, D_hl,
        (1-E_k)*(1-FF_l), E_o, Tp_o, Te_o,
        (1-E_l)*FF_l, Te, De,
        dt=dt, tol=dose_tolerance)
    # Liver
    cb = 1000*J_aorta/CO # mM
    #ca = pkmods.flux_pfcomp(cb, Te, De, dt=dt, solver='interp')
    #Ce = ve*ca/(1-Hct)
    #Ch = dcmri.conc_comp(k_he*ca, Th, t)
    J_aorta = pkmods.flux_pfcomp(J_aorta, Te, De, dt=dt, solver='interp')
    Ne = Tel*FF_l*J_aorta
    Nh = dcmri.conc_comp(E_l*FF_l*J_aorta, Th, t)
    Ce, Ch = Ne/liver_volume, Nh/liver_volume
    if return_conc:
        return t, cb, Ce, Ch
    # Return R
    rp = const.rp(field_strength)
    rh = const.rh(field_strength)
    R1l = R10l + rp*Ce + rh*Ch
    R1b = R10b + rp*cb
    Sb = dcmri.signal_spgress(TR, FA, R1b, S0b)
    Sl = dcmri.signal_spgress(TR, FA, R1l, S0l)
    if not sample:
        return t, Sb, Sl
    Sb = dcmri.sample(t, Sb, xdata, xdata[1]-xdata[0])
    Sl = dcmri.sample(t, Sl, xdata, xdata[1]-xdata[0])
    return np.concatenate([Sb, Sl])


def init(xdata, ydata, pars, vars):

    # Adjust dt if needed
    duration = vars['weight']*vars['dose']/vars['rate']
    if duration != 0:
        if vars['dt'] > duration/5:
            vars['dt'] = duration/5

    n = xdata.size/2
    x = xdata[:n]
    yb = ydata[:n]
    yl = ydata[n:]

    TR, FA = vars['TR'], vars['FA']

    # Estimate BAT and aorta S0 from data
    BAT = x[np.argmax(yb)]
    baseline = x[x <= BAT-20]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    Sref = dcmri.signal_spgress(TR, FA, vars['R10b'], 1)
    S0b = np.mean(yb[:n0]) / Sref

    # Estimate liver S0 from data
    Sref = dcmri.signal_spgress(TR, FA, vars['R10l'], 1)
    S0l = np.mean(yl[:n0]) / Sref

    pars[0] = BAT
    vars['S0b'] = S0b
    vars['S0l'] = S0l
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
    ['E_k', "Kidney extraction fraction", 0.05,"", 0.01, 0.15, 0],
    # Liver
    ['T_g', "Gut mean transit time", 10.0, "sec", 1, 20, 0],
    ['FF_l', "Liver flow fraction", 0.25, "", 0.01, 1, 0],
    ['E_l', "Liver extraction fraction", 0.04, "", 0.00, 0.10, 0],
    ['Te', "Extracellular transit time", 30,"sec", 15, 60.0, 0],
    ['Th', "Hepatocellular transit time", 30*60, "sec", 10*60, 600*60, 0],
]


def export_pars(pars, liver_volume, CO):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters 
    Kve = 1/p.Te
    Ke = Kve/(1-p.E_l)
    Khe = p.E_l*Ke
    lev = p.FF_l*p.CO*p.Te/liver_volume
    Kbh = 1/p.Th
    df.loc['Khe'] = ["Hepatocellular tissue uptake rate", 6000*Khe, 'mL/min/100mL']         
    df.loc['Kbh'] = ["Hepatocellular excretion rate", 6000*Kbh, 'mL/min/100mL']
    df.loc['LBF'] = ["Liver blood flow", 60*p.FF_l*CO, 'mL/min']
    df.loc['LP'] = ["Liver perfusion", 6000*p.FF_l*CO/liver_volume, 'mL/min/100mL']
    df.loc['LEHV'] = ["Liver extra-hepatocellular volume", 100*lev, 'mL/100mL']
    df.loc['CL_l'] = ['Liver blood clearance', 60*p.E_l*p.FF_l*CO, 'mL/min']
    df.loc['E_lb'] = ['Liver whole body extraction fraction', 100*p.E_l*p.FF_l, '%']
    df.loc['k_he'] = ["Hepatocellular tissue uptake rate", 6000*Khe*lev, 'mL/min/100mL'] 
    df.loc['k_bh'] = ["Hepatocellular excretion rate", 6000*Kbh*(1-lev), 'mL/min/100mL'] 
    # Convert to conventional units
    df.loc['FF_l', ['value', 'unit']] = [100*p.FF_l, '%']
    df.loc['E_l', ['value', 'unit']] = [100*p.E_l, '%']
    df.loc['Th', ['value', 'unit']] = [p.Th/60, 'min']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    # This needs a new function plotting aorta and liver
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
