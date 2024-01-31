import numpy as np

def K_flow_1d(dx, u):
    nc = len(u)-1
    # Calculate Kn
    Kn = np.zeros(nc)
    un = u[:-1]
    neg = np.where(un < 0)
    Kn[neg] = -un[neg]/dx
    # Calculate Kp
    Kp = np.zeros(nc)
    up = u[1:]
    pos = np.where(up > 0)
    Kp[pos] = up[pos]/dx     
    return Kp, Kn


def K_diff_1d(dx, D):
    Kn = D[:-1]/dx**2
    Kp = D[1:]/dx**2   
    return Kp, Kn


def K_flowdiff_1d(dx, u, D):
    Ku = K_flow_1d(dx, u)
    Kd = K_diff_1d(dx, D)
    Kp = Ku[0] + Kd[0]
    Kn = Ku[1] + Kd[1]
    return Kp, Kn


def conc_1d1c(t, Jp, Kp, Kn):
    """Concentration in a spatial 1-compartment model in 1D"""
    nt = len(Jp)
    nc = len(Kp)
    K = Kp + Kn
    C = np.zeros((nt,nc))
    for k in range(nt-1):
        dt = t[k+1] - t[k]
        # Initialise at current concentration
        C[k+1,:] = C[k,:]
        # Add influxes at the boundaries:
        C[k+1,0] += dt*(Jp[k+1]+Jp[k])/2
        # Remove outflux to the neigbours:
        C[k+1,:] -= dt*K*C[k,:]
        # Add influx from the neighbours:
        C[k+1,:-1] += dt*Kn[1:]*C[k,1:]
        C[k+1,1:] += dt*Kp[:-1]*C[k,:-1]
    return C


def dt_1d1c(Kp, Kn):
    """maximal time step"""
    K = Kp + Kn
    return np.amin(1/K)

def conc_1d2cfp(t, Jp1, Kp1, Kn1, K02, E21):
    """Concentration in a spatial 2-compartment filtration model in 1D (positive influx only)"""
    nt = len(Jp1)
    nc = len(Kp1)
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = (Kp1 + Kn1)*E21/(1-E21)
    K1 = Kp1 + Kn1 + K21
    K2 = K02
    C1 = np.zeros((nt,nc))
    C2 = np.zeros((nt,nc))
    for k in range(nt-1):
        dt = t[k+1]-t[k]
        # Initialise at current concentration
        C1[k+1,:] = C1[k,:]
        C2[k+1,:] = C2[k,:]
        # Add influxes at the boundaries:
        C1[k+1,0] += dt*Jp1[k]
        # Remove outflux to the neigbours:
        C1[k+1,:] -= dt*K1*C1[k,:]
        C2[k+1,:] -= dt*K2*C2[k,:]
        # Add influx from the neighbours:
        C1[k+1,:-1] += dt*Kn1[1:]*C1[k,1:]
        C1[k+1,1:] += dt*Kp1[:-1]*C1[k,:-1]
        # Add influx at the same location
        C2[k+1,:] += dt*K21*C1[k,:]
    return C1, C2


def dt_1d2cfp(Kp1, Kn1, K02, E21):
    """maximal time step"""
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = (Kp1 + Kn1)*E21/(1-E21)
    K1 = Kp1 + Kn1 + K21
    K2 = K02
    # dt*K1<1
    # dt*K2<1
    return np.amin([np.amin(1/K1), np.amin(1/K2)])


def conc_1d3cf(t, Jp1, Kp1, K32, E21, Kn3):
    """Concentration in a spatial 3-compartment filtration model in 1D (positive influx only)"""
    nt = len(Jp1)
    nc = len(Kp1)
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = Kp1*E21/(1-E21)
    K1 = Kp1 + K21
    K2 = K32
    K3 = Kn3
    C1 = np.zeros((nt,nc))
    C2 = np.zeros((nt,nc))
    C3 = np.zeros((nt,nc))
    for k in range(nt-1):
        dt = t[k+1]-t[k]
        # Initialise at current concentration
        C1[k+1,:] = C1[k,:]
        C2[k+1,:] = C2[k,:]
        C3[k+1,:] = C3[k,:]
        # Add influxes at the boundaries:
        C1[k+1,0] += dt*Jp1[k]
        # Remove outflux to the neigbours:
        C1[k+1,:] -= dt*K1*C1[k,:]
        C2[k+1,:] -= dt*K2*C2[k,:]
        C3[k+1,:] -= dt*K3*C3[k,:]
        # Add influx from the neighbours:
        C1[k+1,1:] += dt*Kp1[:-1]*C1[k,:-1]
        C3[k+1,:-1] += dt*Kn3[1:]*C3[k,1:]
        # Add influx at the same location
        C2[k+1,:] += dt*K21*C1[k,:]
        C3[k+1,:] += dt*K32*C2[k,:]
    return C1, C2, C3


def dt_1d3cf(Kp1, K32, E21, Kn3):
    """maximal time step"""
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = Kp1*E21/(1-E21)
    K1 = Kp1 + K21
    K2 = K32
    K3 = Kn3
    # dt*K1<1
    # dt*K2<1
    return np.amin([np.amin(1/K1), np.amin(1/K2), np.amin(1/K3)])


def conc_1d2cxp(t, Jp1, Kp1, Kn1, K12, E21):
    """Concentration in a spatial 2-compartment exchange model in 1D (positive influx only)"""
    nt = len(Jp1)
    nc = len(Kp1)
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = (Kp1 + Kn1)*E21/(1-E21)
    K1 = Kp1 + Kn1 + K21
    K2 = K12
    C1 = np.zeros((nt,nc))
    C2 = np.zeros((nt,nc))
    for k in range(nt-1):
        dt = t[k+1]-t[k]
        # Initialise at current concentration
        C1[k+1,:] = C1[k,:]
        C2[k+1,:] = C2[k,:]
        # Add influxes at the boundaries:
        C1[k+1,0] += dt*Jp1[k]
        # Remove outflux to the neigbours:
        C1[k+1,:] -= dt*K1*C1[k,:]
        C2[k+1,:] -= dt*K2*C2[k,:]
        # Add influx from the neighbours:
        C1[k+1,:-1] += dt*Kn1[1:]*C1[k,1:]
        C1[k+1,1:] += dt*Kp1[:-1]*C1[k,:-1]
        # Add influx at the same location
        C2[k+1,:] += dt*K21*C1[k,:]
        C1[k+1,:] += dt*K12*C2[k,:]
    return C1, C2


def dt_1d2cxp(Kp1, Kn1, K12, E21):
    """maximal time step"""
    # K21 = E21*(Kp1 + Kn1 + K21)
    # K21(1-E21) = E21*(Kp1 + Kn1)
    K21 = (Kp1 + Kn1)*E21/(1-E21)
    K1 = Kp1 + Kn1 + K21
    K2 = K12
    # dt*K1<1
    # dt*K2<1
    return np.amin([np.amin(1/K1), np.amin(1/K2)])
