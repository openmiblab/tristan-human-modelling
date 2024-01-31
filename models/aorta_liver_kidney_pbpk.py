import numpy as np

import dcmri
import dcmods.pkmods as pkmods
import dcmods.const as const
import dcmods.tools as tools
import dcmods.fig as fig


def model(xdata,
        BAT, CO, T_hl, D_hl, E_o, Tp_o, Te_o,
        T_g, FF_l, E_l, Te, Th,
        Ta, FF_k, F_k, Tp, dr, ut0, ut1, ut2, r0, r1, r2,
        R10l = 1,
        R10b = 1,
        R10k = 1,
        S0l = 1,
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
        liver_volume = None,
        return_conc = False,
        sample = True,
        ):
    t = np.arange(0, max(xdata)+xdata[1]+dt, dt)
    # Aorta
    Ji = dcmri.inj_step(t, # mmol/sec
        weight, conc, dose, rate, BAT)
    J_aorta = pkmods.prop_body2_liver_kidneys(t, Ji,
        T_hl, D_hl,
        E_o, Tp_o, Te_o,
        T_g, FF_l, E_l, 1/Te, 
        FF_k, E_k, 1/Tp,
        tol=dose_tolerance)
    cb = 1000*J_aorta/CO # mM
    # Liver
    J_liver = FF_l*J_aorta
    J_liver = dcmri.flux_comp(J_liver, T_g, t)
    Ne, Nh = pkmods.conc_liver(t, J_liver, 1/Te, 0, 1/Th, E_l)
    Ce_l = 1000*Ne/liver_volume # mM
    Ch_l = 1000*Nh/liver_volume # mM
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
        return t, cb, Ce_l, Ch_l, Cp, Ct
    # Return R
    rp = const.rp(field_strength)
    rh = const.rh(field_strength)
    R1l = R10l + rp*Ce_l + rh*Ch_l
    R1k = R10k + rp*Cp + rp*Ct
    R1b = R10b + rp*cb
    Sb = dcmri.signal_spgress(TR, FA, R1b, S0b)
    Sl = dcmri.signal_spgress(TR, FA, R1l, S0l)
    Sk = dcmri.signal_spgress(TR, FA, R1k, S0k)
    if not sample:
        return t, Sb, Sl, Sk
    Sb = dcmri.sample(t, Sb, xdata, xdata[1]-xdata[0])
    Sl = dcmri.sample(t, Sl, xdata, xdata[1]-xdata[0])
    Sk = dcmri.sample(t, Sk, xdata, xdata[1]-xdata[0])
    return np.concatenate([Sb, Sl, Sk])


def init(xdata, ydata, pars, vars):

    # Adjust dt if needed
    duration = vars['weight']*vars['dose']/vars['rate']
    if duration != 0:
        if vars['dt'] > duration/5:
            vars['dt'] = duration/5

    n = xdata.size/3
    xb, yb = xdata[:n], ydata[:n]
    xl, yl = xdata[n:2*n], ydata[n:2*n]
    xk, yk = xdata[2*n:], ydata[2*n:]

    TR, FA = vars['TR'], vars['FA']

    # Estimate BAT from data
    BAT = xb[np.argmax(yb)]
    baseline = xb[xb <= BAT-20]
    n0 = baseline.size
    if n0 == 0: 
        n0 = 1
    pars[0] = BAT

    # Estimate S0 from data
    Sref = dcmri.signal_spgress(TR, FA, vars['R10b'], 1)
    vars['S0b'] = np.mean(yb[:n0]) / Sref
    Sref = dcmri.signal_spgress(TR, FA, vars['R10l'], 1)
    vars['S0l'] = np.mean(yl[:n0]) / Sref
    Sref = dcmri.signal_spgress(TR, FA, vars['R10k'], 1)
    vars['S0k'] = np.mean(yk[:n0]) / Sref

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
    # Liver
    ['T_g', "Gut mean transit time", 10.0, "sec", 1, 20, 0],
    ['FF_l', "Liver flow fraction", 0.25, "", 0.01, 1, 0],
    ['E_l', "Liver extraction fraction", 0.04, "", 0.00, 0.10, 0],
    ['Te', "Extracellular transit time", 30,"sec", 15, 60.0, 0],
    ['Th', "Hepatocellular transit time", 30*60, "sec", 10*60, 600*60, 0],
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


def export_pars(pars, kidney_volume, liver_volume, CO, Hct):
    df = tools.df_export(PARS, pars)
    p = df.value
    # Add derived parameters - kidney
    E_k = np.divide(p.F_k, 1+p.F_k)
    df.loc['E_k'] = ["Kidney extraction fraction", 100*E_k, '%']
    df.loc['GFR'] = ["Glomerular Filtration Rate", 60*p.F_k*p.FF_k*CO*(1-Hct), 'mL/min']  
    df.loc['RBF'] = ["Renal blood flow", 60*p.FF_k*CO, 'mL/min']
    df.loc['RP'] = ["Renal perfusion", 6000*p.FF_k*CO/kidney_volume, 'mL/min/100mL']
    df.loc['RBV'] = ["Renal blood volume", 100*p.FF_k*CO*p.Tv/kidney_volume, 'mL/100mL']
    # Convert to conventional units - kidney
    df.loc['E_o', ['value', 'unit']] = [100*p.E_o, '%']
    df.loc['E_l', ['value', 'unit']] = [100*p.E_l, '%']
    df.loc['FF_k', ['value', 'unit']] = [100*p.FF_k, '%']
    df.loc['F_k', ['value', 'unit']] = [100*p.F_k, '%']
    # Add derived parameters - liver
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
    # Convert to conventional units - liver
    df.loc['FF_l', ['value', 'unit']] = [100*p.FF_l, '%']
    df.loc['E_l', ['value', 'unit']] = [100*p.E_l, '%']
    df.loc['Th', ['value', 'unit']] = [p.Th/60, 'min']
    return df


def plot_fit(xdata, ydata, pars, BAT, vars={}, xvalid=None, show=True, save=False, path=None, prefix=''):
    clrs = ['cornflowerblue','darkblue'] 
    label = ['Plasma', 'Tubuli']
    # This needs a new function plotting aorta/liver/kidney
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win', xlim=[BAT-20, BAT+40], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
    fig.tissue_2cm(model, xdata, ydata, pars, vars=vars, xvalid=xvalid, win='win_', xlim=[BAT-20, BAT+160], show=show, save=save, path=path, prefix=prefix, color=clrs, label=label)
