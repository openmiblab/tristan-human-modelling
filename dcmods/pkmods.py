import numpy as np
from scipy.integrate import trapezoid
import dcmri
import dcmods.pkcont as pkcont

# Blocks


def flux_pfcomp(J, T, D, t=None, dt=1.0, solver='conv'):
    if D == 0:
        return dcmri.flux_plug(J, T, t=t, dt=dt, solver=solver)
    if D == 1:
        return dcmri.flux_comp(J, T, t=t, dt=dt)
    Tc = D*T
    Tp = (1-D)*T
    J = dcmri.flux_comp(J, Tc, t=t, dt=dt)
    J = dcmri.flux_plug(J, Tp, t=t, dt=dt, solver=solver)
    return J



# Whole body


def flux_wb_cc(J_lungs,
        T_lh, 
        E_organs, Tp_organs, Te_organs,
        E_extr,
        t=None, dt=1.0, tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, x=t, dx=dt)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_comp(J_lungs, T_lh, t=t, dt=dt)
        # Propagate through the other organs
        J_aorta = np.stack((J_aorta, np.zeros(J_aorta.size)))
        J_lungs = (1-E_extr)*dcmri.flux_2comp(J_aorta, [Tp_organs, Te_organs], E_o, t=t, dt=dt)[0,0,:]
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, x=t, dx=dt)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_comp(J_lungs_total, T_lh, t=t, dt=dt)
    # Return total flux into aorta
    return J_aorta_total


def flux_wb_pcc(J_lungs,
        T_lh, D_lh,
        E_organs, Tp_organs, Te_organs,
        E_extr,
        t=None, dt=1.0, tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, x=t, dx=dt)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = flux_pfcomp(J_lungs, T_lh, D_lh, t=t, dt=dt, solver='interp')
        # Propagate through the other organs
        J_aorta = np.stack((J_aorta, np.zeros(J_aorta.size)))
        J_lungs = (1-E_extr)*dcmri.flux_2comp(J_aorta, [Tp_organs, Te_organs], E_o, t=t, dt=dt)[0,0,:]
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, x=t, dx=dt)
    # Propagate total flux through lungs
    J_aorta_total = flux_pfcomp(J_lungs_total, T_lh, D_lh, t=t, dt=dt, solver='interp')
    # Return total flux into aorta
    return J_aorta_total


def flux_wb_chc(J_lungs,
        T_lh, D_lh,
        E_organs, Tp_organs, Te_organs,
        E_extr,
        t=None, dt=1.0, tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, x=t, dx=dt)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_chain(J_lungs, T_lh, D_lh, t=t, dt=dt, solver='step')
        # Propagate through the other organs
        J_aorta = np.stack((J_aorta, np.zeros(J_aorta.size)))
        J_lungs = (1-E_extr)*dcmri.flux_2comp(J_aorta, [Tp_organs, Te_organs], E_o, t=t, dt=dt)[0,0,:]
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, x=t, dx=dt)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_chain(J_lungs_total, T_lh, D_lh, t=t, dt=dt, solver='step')
    # Return total flux into aorta
    return J_aorta_total


def flux_wb_liver(J_lungs,
        T_lh, D_lh,
        R_organs, E_organs, Tp_organs, Te_organs,
        R_liver, Te_liver, De_liver,
        t=None, dt=1.0, tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, x=t, dx=dt)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = flux_pfcomp(J_lungs, T_lh, D_lh, t=t, dt=dt)
        # Propagate through liver and other organs
        # R_liver = (1-E_liver)*FF_liver
        # R_organs = (1-E_kidneys)*(1-FF_liver)
        J_liver = R_liver*flux_pfcomp(J_aorta, Te_liver, De_liver, t=t, dt=dt)
        J_aorta = np.stack((J_aorta, np.zeros(J_aorta.size)))
        J_organs = R_organs*dcmri.flux_2comp(J_aorta, [Tp_organs, Te_organs], E_o, t=t, dt=dt)[0,0,:]
        # Add up outfluxes from liver and other organs
        J_lungs = J_liver + J_organs
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, x=t, dx=dt)
    # Propagate total flux through lungs
    J_aorta_total = flux_pfcomp(J_lungs_total, T_lh, D_lh, t=t, dt=dt)
    # Return total flux into aorta
    return J_aorta_total


def prop_body_kidneys(t, J_lungs,
        T_lh, 
        E_organs, Tp_organs, Te_organs,
        E_liver,
        FF_kidneys, E_kidneys, Ke_kidneys, 
        tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, t)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_comp(J_lungs, T_lh, t)
        # Split into liver and other organs
        J_kidneys = FF_kidneys*J_aorta
        J_organs = (1-FF_kidneys)*J_aorta
        # Propagate through liver and other organs
        J_kidneys = (1-E_kidneys)*dcmri.flux_nscomp(J_kidneys, 1/Ke_kidneys, t)
        J_organs = (1-E_liver)*dcmri.flux_2comp(J_organs, [Tp_organs, Te_organs], E_o, t)[0,0,:]
        # Add up outfluxes from liver and other organs
        J_lungs = J_kidneys + J_organs
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, t)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_comp(J_lungs_total, T_lh, t)
    # Return total flux into aorta
    return J_aorta_total

def prop_body2_kidneys(t, J_lungs,
        T_l, T_h,
        E_organs, Tp_organs, Te_organs,
        E_liver,
        FF_kidneys, E_kidneys, Ke_kidneys, 
        T_venous,
        tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, t)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_comp(J_lungs, T_l, t)
        J_aorta = dcmri.flux_comp(J_aorta, T_h, t)
        # Split into liver and other organs
        J_kidneys = FF_kidneys*J_aorta
        J_organs = (1-FF_kidneys)*J_aorta
        # Propagate through liver and other organs
        J_kidneys = (1-E_kidneys)*dcmri.flux_nscomp(J_kidneys, 1/Ke_kidneys, t)
        J_organs = (1-E_liver)*dcmri.flux_2comp(J_organs, [Tp_organs, Te_organs], E_o, t)[0,0,:]
        # Add up outfluxes from liver and other organs
        J_lungs = J_kidneys + J_organs
        # Propgate through venous return
        J_lungs = dcmri.flux_plug(J_lungs, T_venous, t)
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, t)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_comp(J_lungs_total, T_l, t)
    J_aorta_total = dcmri.flux_comp(J_aorta_total, T_h, t)
    # Return total flux into aorta
    return J_aorta_total

def prop_body3_kidneys(t, J_lungs,
        T_lh, D_lh,
        E_organs, Tp_organs, Te_organs,
        E_liver,
        FF_kidneys, E_kidneys, Ke_kidneys, Ta_kidneys,
        tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, t)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_chain(J_lungs, T_lh, D_lh, dt=t[1]-t[0])
        # Split into liver and other organs
        J_kidneys = FF_kidneys*J_aorta
        J_organs = (1-FF_kidneys)*J_aorta
        # Propagate through liver and other organs
        J_kidneys = dcmri.flux_plug(J_kidneys, Ta_kidneys, t)
        J_kidneys = (1-E_kidneys)*dcmri.flux_nscomp(J_kidneys, 1/Ke_kidneys, t)
        J_organs = (1-E_liver)*dcmri.flux_2comp(J_organs, [Tp_organs, Te_organs], E_o, t)[0,0,:]
        # Add up outfluxes from liver and other organs
        J_lungs = J_kidneys + J_organs
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, t)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_chain(J_lungs_total, T_lh, D_lh, dt=t[1]-t[0])
    # Return total flux into aorta
    return J_aorta_total


def prop_body_liver_kidneys(t, J_lungs,
        T_lh, 
        E_organs, Tp_organs, Te_organs,
        T_gut, FF_liver, E_liver, Ke_liver, 
        FF_kidneys, E_kidneys, Kp_kidneys,
        tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, t)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_comp(J_lungs, T_lh, t)
        # Split into liver, kidneys and other organs
        J_liver = FF_liver*J_aorta
        J_kidneys = FF_kidneys*J_aorta
        J_organs = (1-FF_liver-FF_kidneys)*J_aorta
        # Propagate through liver, kidneys and other organs
        J_liver = dcmri.flux_comp(J_liver, T_gut, t)
        J_liver = (1-E_liver)*dcmri.flux_nscomp(J_liver, 1/Ke_liver, t)
        J_kidneys = (1-E_kidneys)*dcmri.flux_nscomp(J_kidneys, Kp_kidneys, t)
        J_organs = dcmri.flux_2comp(J_organs, [Tp_organs, Te_organs], E_o, t)[0,0,:]
        # Add up outfluxes from liver, kidneys and other organs
        J_lungs = J_liver + J_kidneys + J_organs
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, t)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_comp(J_lungs_total, T_lh, t)
    # Return total flux into aorta
    return J_aorta_total


def prop_body2_liver_kidneys(t, J_lungs,
        T_lh, D_lh,
        E_organs, Tp_organs, Te_organs,
        T_gut, FF_liver, E_liver, Ke_liver, 
        FF_kidneys, E_kidneys, Kp_kidneys,
        tol = 0.001):
    E_o = [[1-E_organs,1],[E_organs,0]]
    dose0 = np.trapz(J_lungs, t)
    dose = dose0
    min_dose = tol*dose0
    J_lungs_total = J_lungs
    while dose > min_dose:
        # Propagate through the lungs and heart
        J_aorta = dcmri.flux_chain(J_lungs, T_lh, D_lh, dt=t[1]-t[0])
        # Split into liver, kidneys and other organs
        J_liver = FF_liver*J_aorta
        J_kidneys = FF_kidneys*J_aorta
        J_organs = (1-FF_liver-FF_kidneys)*J_aorta
        # Propagate through liver, kidneys and other organs
        J_liver = dcmri.flux_comp(J_liver, T_gut, t)
        J_liver = (1-E_liver)*dcmri.flux_nscomp(J_liver, 1/Ke_liver, t)
        J_kidneys = (1-E_kidneys)*dcmri.flux_nscomp(J_kidneys, Kp_kidneys, t)
        J_organs = dcmri.flux_2comp(J_organs, [Tp_organs, Te_organs], E_o, t)[0,0,:]
        # Add up outfluxes from liver, kidneys and other organs
        J_lungs = J_liver + J_kidneys + J_organs
        # Add to the total flux into the lungs
        J_lungs_total += J_lungs
        # Get residual dose
        dose = np.trapz(J_lungs, t)
    # Propagate total flux through lungs
    J_aorta_total = dcmri.flux_chain(J_lungs_total, T_lh, D_lh, t[1]-t[0])
    # Return total flux into aorta
    return J_aorta_total


# Liver


def conc_liver(
            # Liver model 
            t, # time points (sec)
            J, # influx per unit volume (mmol/sec/mL tissue)
            vel, # velocity (1/sec)
            dif, # diffusion rate (1/sec)
            Kbh, # Biliary excretion rate (sec)
            El, # Extraction fraction
            n = 20, # nr of numerical voxels
            ):
    # dimensionless locations of numerical voxel boundaries
    x = np.linspace(0, 1, n+1)
    dx = x[1] 
    # interpolate parameters on tube 
    vel_x = dcmri.interp(x, vel, pos=True)
    dif_x = dcmri.interp(x, dif, pos=True)
    Kbh_x = dcmri.interp(x, Kbh, pos=True)
    El_x = dcmri.interp(x, El, pos=True)
    Kbh_x = (Kbh_x[1:]+Kbh_x[:-1])/2
    El_x = (El_x[1:]+El_x[:-1])/2
    # Compartmental rate constants
    Kp, Kn = pkcont.K_flowdiff_1d(dx, vel_x, dif_x) # 1/sec
    # High resolution time points
    #dth = 0.9/(2*dif/dx**2 + np.amax(vel_x)/dx)
    dth = 0.9*pkcont.dt_1d2cfp(Kp, Kn, Kbh_x, El_x)
    tacq = np.amax(t)
    if tacq/dth > 1e6:
        print('dth', dth, tacq/dth)
        print('vel', vel)
        print('vel_x', dcmri.interp(x, vel, pos=True))
        print('dif', dif)
        print('dif_x', dif_x)
        print('Ti', Kbh)
        print('Ti_x', Kbh_x)
        print('Ei', El)
        print('Ei_x', El_x)
        print('Kp', Kp)
        print('Kn', Kn)
        msg = 'Rate constants are too high for numerical computation. Consider reducing their upper bound.'
        raise ValueError(msg)
    nth = 1 + np.ceil(tacq/dth).astype(np.int32)
    th = np.linspace(0, tacq, nth)
    # Upsample influx
    Jh = np.interp(th, t, J/dx)
    # Calculate concentrations
    Ce, Ch = pkcont.conc_1d2cfp(th, Jh, Kp, Kn, Kbh_x, El_x)
    # Sum concentrations over all voxels
    Ce = dx*np.sum(Ce, axis=1)
    Ch = dx*np.sum(Ch, axis=1)
    # Downsample concentrations to measured time resolution.
    Ce = np.interp(t, th, Ce)
    Ch = np.interp(t, th, Ch)
    return Ce, Ch # mmol/mL tissue



# Kidney


def conc_neph(
            # Nephron model with linearly varying reabsorption
            t, # time points (sec)
            J, # influx per unit volume (mmol/sec/mL tissue)
            vel, # velocity (1/sec)
            dif, # diffusion rate (1/sec)
            reabs, # reabsorption rate
            n = 20, # nr of numerical voxels
            ):
    # dimensionless locations of numerical voxel boundaries
    x = np.linspace(0, 1, n+1)
    dx = x[1] 
    # reabsorption rate everywhere
    reabs_x = dcmri.interp(x, reabs, pos=True)
    # velocity everywhere
    vel_x = dcmri.interp(x, vel, pos=True)
    vel_x = vel_x*np.exp(-trapezoid(reabs_x, dx=dx)) # reabs may be unnecessary
    # diffusion rate everywhere 
    dif_x = dcmri.interp(x, dif, pos=True)
    # Compartmental rate constants
    Kp, Kn = pkcont.K_flowdiff_1d(dx, vel_x, dif_x) # 1/sec
    # High resolution time points
    #dth = 0.9/(2*dif/dx**2 + np.amax(vel_x)/dx)
    dth = 0.9 * pkcont.dt_1d1c(Kp, Kn)
    tacq = np.amax(t)
    if tacq/dth > 1e6:
        print('dth', dth, tacq/dth)
        print('vel', vel)
        print('vel_x', dcmri.interp(x, vel, pos=True))
        print('dif', dif)
        print('dif_x', dif_x)
        print('reabs', reabs)
        print('reabs_x', reabs_x)
        print('Kp', Kp)
        print('Kn', Kn)
        msg = 'Tubular velocity is too high. Reduce upper bound on the velocity.'
        raise ValueError(msg)
    nth = 1 + np.ceil(tacq/dth).astype(np.int32)
    th = np.linspace(0, tacq, nth)
    # Upsample influx
    Jh = np.interp(th, t, J/dx)
    # Calculate concentrations
    Ch = pkcont.conc_1d1c(th, Jh, Kp, Kn)
    # Sum concentrations over all voxels
    Cth = dx*np.sum(Ch, axis=1)
    # Downsample concentrations to measured time resolution
    Ct = np.interp(t, th, Cth)
    return Ct # mmol/mL tissue


def conc_nephi(
            # Nephron model with linearly varying reabsorption
            t, # time points (sec)
            J, # influx per unit volume (mmol/sec/mL tissue)
            vel, # velocity (1/sec)
            dif, # diffusion rate (1/sec)
            reabs, # reabsorption rate
            Ti, # Interstitial MTT (sec)
            Ei, # Interstitial extraction fraction
            n = 20, # nr of numerical voxels
            ):
    # dimensionless locations of numerical voxel boundaries
    x = np.linspace(0, 1, n+1)
    dx = x[1] 
    # interpolate parameters on tube 
    vel_x = dcmri.interp(x, vel, pos=True)
    dif_x = dcmri.interp(x, dif, pos=True)
    reabs_x = dcmri.interp(x, reabs, pos=True)
    Ti_x = dcmri.interp(x, Ti, pos=True)
    Ei_x = dcmri.interp(x, Ei, pos=True)
    Ti_x = (Ti_x[1:]+Ti_x[:-1])/2
    Ei_x = (Ei_x[1:]+Ei_x[:-1])/2
    # reabsorption correction
    vel_x = vel_x*np.exp(-trapezoid(reabs_x, dx=dx))
    # Compartmental rate constants
    Kp, Kn = pkcont.K_flowdiff_1d(dx, vel_x, dif_x) # 1/sec
    # High resolution time points
    #dth = 0.9/(2*dif/dx**2 + np.amax(vel_x)/dx)
    dth = 0.9*pkcont.dt_1d2cxp(Kp, Kn, 1/Ti_x, Ei_x)
    tacq = np.amax(t)
    if tacq/dth > 1e6:
        print('dth', dth, tacq/dth)
        print('vel', vel)
        print('vel_x', dcmri.interp(x, vel, pos=True))
        print('dif', dif)
        print('dif_x', dif_x)
        print('reabs', reabs)
        print('reabs_x', reabs_x)
        print('Ti', Ti)
        print('Ti_x', Ti_x)
        print('Ei', Ei)
        print('Ei_x', Ei_x)
        print('Kp', Kp)
        print('Kn', Kn)
        msg = 'Rate constants are too high for numerical computation. Consider reducing their upper bound.'
        raise ValueError(msg)
    nth = 1 + np.ceil(tacq/dth).astype(np.int32)
    th = np.linspace(0, tacq, nth)
    # Upsample influx
    Jh = np.interp(th, t, J/dx)
    # Calculate concentrations
    # Ch = conc_1d1c(th, Jh, Kp, Kn)
    Cht, Chi = pkcont.conc_1d2cxp(th, Jh, Kp, Kn, 1/Ti_x, Ei_x)
    Ch = Cht + Chi
    # Sum concentrations over all voxels
    Cth = dx*np.sum(Ch, axis=1)
    # Downsample concentrations to measured time resolution.
    Ct = np.interp(t, th, Cth)
    return Ct # mmol/mL tissue
