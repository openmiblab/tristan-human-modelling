import warnings
import os
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit as scipy_curve_fit

def read_csv(file):
    try:
        df = pd.read_csv(file, index_col='symbol')   
    except:
        msg = 'Cannot read model parameters from file ' + file
        msg += '\nPlease check if the file exists and is not open in another program.'
        raise RuntimeError(msg)
    return df

def df_export(PARS, pars):
    df = df_pars(PARS, pars)
    df = df.drop(['initial value','lower bound','upper bound','fix'], axis=1)
    return df

def df_pars(PARS, pars):
    cols = ['symbol', "name", "initial value", "unit", "lower bound", "upper bound", "fix"]
    p = pd.DataFrame(PARS, columns=cols)
    p.set_index('symbol', inplace=True)
    p['value'] = p['initial value']
    cols = cols[1:]
    cols.insert(2, 'value')
    p = p[cols]
    p.loc[:,'value'] = pars
    return p

def to_csv(PARS, vals, file):
    p = df_pars(PARS, vals)
    path = os.path.dirname(file)
    if not os.path.isdir(path):
        os.makedirs(path)
    try:
        p.to_csv(file)
    except:
        print("Can't write to file ", file)
        print("Please close the file before saving data.") 


def unpack_pars(p):
    p0 = [par[2] for par in p]
    b0 = [par[4] for par in p]
    b1 = [par[5] for par in p]
    fix = [par[6] for par in p]
    return p0, (b0, b1), fix


def save_path(path):
    if path is None:
        path = os.path.dirname(__file__)
        path = os.path.join(path, 'results')
        if not os.path.isdir(path):
            os.mkdir(path)
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def xfit(xdata, xrange=None, xvalid=None):
    # Identify x- and y-values that are fitted
    if xrange is None:
        xrange = [xdata.min(), xdata.max()]
    if xvalid is None:
        xvalid = np.ones(len(xdata), dtype=np.int32)
    elif not isinstance(xvalid, np.ndarray):
        xvalid = np.array(xvalid, dtype=np.int32)
    xfit = xvalid * (xrange[0]<=xdata) * (xdata<=xrange[1])
    return xfit==1


def loss(f, xdata, ydata, p0, vars={}, xrange=None, xvalid=None):
    # Predict data at all xvalues
    y = f(xdata, *p0, **vars)
    # Identify x- and y-values that are fitted
    xf = xfit(xdata, xrange, xvalid)
    # Calclulate the loss at the fitted values
    loss = np.linalg.norm(y[xf] - ydata[xf])/np.linalg.norm(ydata[xf])
    return 100*loss


# f(xdata, *pars, **kwargs)
# wrapper for curve_fit except that p0 is REQUIRED
def curve_fit(f, xdata, ydata, p0, vars={}, pfix=None, xrange=None, xvalid=None,  sigma=None, absolute_sigma=False, check_finite=None, bounds=(-np.inf, np.inf), method=None, jac=None, full_output=False, nan_policy=None, **kwargs):
    
    # Convert array_like to array
    if not isinstance(p0, np.ndarray):
        p0 = np.array(p0).astype(np.float64)

    # Identify parameters that are fitted
    if pfix is None:
        pfit = np.full(len(p0), True)
    else:
        if not isinstance(pfix, np.ndarray):
            pfix = np.array(pfix)
        if np.sum(pfix) == pfix.size:
            msg = 'All parameters are fixed -- returning initial values!!'
            warnings.warn(msg)
            if full_output:
                return p0, None, None
            else:
                return p0, None
        pfit = pfix==0
        
    # Select bounds for the parameters that are fitted
    b0 = bounds[0]
    if not np.isscalar(b0):
        if not isinstance(b0, np.ndarray):
            b0 = np.array(b0).astype(np.float64)
        b0 = b0[pfit]  
    b1 = bounds[1]
    if not np.isscalar(b1):
        if not isinstance(b1, np.ndarray):
            b1 = np.array(b1).astype(np.float64)
        b1 = b1[pfit]

    # Define a fit function for the parameters that are fitted 
    p0c = p0.copy()
    def fit_func(x,*p):
        p0c[pfit] = p
        return f(x, *p0c, **vars)      

    # Identify x- and y-values that are fitted
    xf = xfit(xdata, xrange, xvalid)

    # Apply curve_fit to restricted data
    try:
        p = scipy_curve_fit(
            fit_func, xdata[xf], ydata[xf], 
            p0 = p0[pfit], 
            sigma = sigma,
            absolute_sigma = absolute_sigma,
            check_finite = check_finite,
            bounds = (b0,b1),
            method = method,
            jac = jac, 
            full_output = full_output,
            nan_policy = nan_policy,
            **kwargs)
    except RuntimeError as e:
        msg = 'Runtime error in curve_fit -- \n'
        msg += str(e) + ' Returning initial values.'
        warnings.warn(msg)
        p = p0
    
    # Return parameter array in original length
    p0c[pfit] = p[0]
    if full_output:
        return (p0c,) + (p[1],) + p[2:]
    else:
        return p0c, p[1]


def test_curve_fit():
    def third_order(x, a, b, c, d, add=2):
        return a + b*x + c*x**2 + (d+add)*x**3
    nx = 5
    x = np.linspace(0,1,nx)
    y = third_order(x, 2, 3, 4, 5, add=5)
    p0 = [1,1,1,1]
    p, pcov = curve_fit(third_order, x, y, 
        p0, 
        kwargs = {'add':5},
        pfix = [0,0,0,0],
        xrange = [0,1.0], # with nx=5 and xrange = [0.3,1.0] this runs into RuntimeError
        xvalid = [1,1,0,0,1],
        bounds = (
            [0,0,0,0],
            [6,6,6,6],
        ), 
    )
    print(p)
    

if __name__ == '__main__':

    test_curve_fit()