import os
import numpy as np

import matplotlib.pyplot as plt
import dcmods.tools as tools


def conc_tissue(ax, xlim, f, xdata, pars, vars={}, color='black', label='Conc'):
    t, cb = f(xdata, *pars, return_conc=True, **vars)
    if xlim is None:
        xlim = [t[0],t[-1]]
    ax.set(xlabel='Time (min)', ylabel='Concentration (mM)', xlim=np.array(xlim)/60)
    ax.plot(t/60, 0*t, color='gray')
    ax.plot(t/60, cb, linestyle='-', color=color, linewidth=3.0, label=label)
    ax.legend()


def conc_2cm(ax, xlim, f, xdata, pars, vars={}, color='black', label=['Extracellular','Hepatocyte']):
    t, Ce, Ch = f(xdata, *pars, return_conc=True, **vars)
    if xlim is None:
        xlim = [t[0],t[-1]]
    t = t/60
    ax.set(xlabel='Time (min)', ylabel='Tissue concentration (mM)', xlim=np.array(xlim)/60)
    ax.plot(t, 0*t, color='gray')
    ax.plot(t, Ce, linestyle='-.', color=color, linewidth=2.0, label=label[0])
    ax.plot(t, Ch, linestyle='--', color=color, linewidth=2.0, label=label[1])
    ax.plot(t, Ce+Ch, linestyle='-', color=color, linewidth=3.0, label='Tissue')
    ax.legend()


def data(ax, xlim, f, xdata, ydata, pars, vars={}, xrange=None, xvalid=None, color=['black', 'black'], xcheck=None, ycheck=None):
    xf = tools.xfit(xdata, xrange, xvalid)
    xn = np.logical_not(xf)
    tacq = xdata[1]-xdata[0]
    t, sig = f(xdata, *pars, sample=False, **vars)
    if xlim is None:
        xlim = [t[0],t[-1]]
    ax.set(xlabel='Time (min)', ylabel='MR Signal (a.u.)', xlim=np.array(xlim)/60)
    ax.plot((xdata[xn]+tacq/2)/60, ydata[xn], marker='o', color='gray', label='ignored data', linestyle = 'None')
    ax.plot((xdata[xf]+tacq/2)/60, ydata[xf], marker='o', color=color[0], label='fitted data', linestyle = 'None')
    ax.plot(t/60, sig, linestyle='-', color=color[1], linewidth=3.0, label='fit' )
    #ax.plot(xdata[-ncal:]/60, ydata[-ncal:], color='black', marker='D', linestyle='None', label='MOLLI')
    if ycheck is not None:
        ax.plot(xcheck/60, ycheck, color='black', marker='D', linestyle='None', label='Checkpoints')
    ax.legend()


def tissue(f, xdata, ydata, pars, vars={}, xvalid=None, win='all', xlim=None, show=True, save=False, path=None, prefix='', color=['black', 'black'], label='Conc', xcheck=None, ycheck=None):
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,4))
    fig.subplots_adjust(wspace=0.3)
    data(ax1, xlim, f, xdata, ydata, pars, vars=vars, xvalid=xvalid, color=color, xcheck=xcheck, ycheck=ycheck)
    conc_tissue(ax2, xlim, f, xdata, pars, vars=vars, color=color[1], label=label)
    if save:   
        path = tools.save_path(path)      
        plt.savefig(fname=os.path.join(path, prefix + '_' + win + '.png'))
    if show:
        plt.show()
    else:
        plt.close()


def tissue_2cm(f, xdata, ydata, pars, vars={}, xvalid=None, win='all', xlim=None, show=True, save=False, path=None, prefix='', color=['black', 'black'], label=['Extracellular','Hepatocyte'], xcheck=None, ycheck=None):
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(10,4))
    fig.subplots_adjust(wspace=0.3)
    data(ax1, xlim, f, xdata, ydata, pars, vars=vars, xvalid=xvalid, color=color, xcheck=xcheck, ycheck=ycheck)
    conc_2cm(ax2, xlim, f, xdata, pars, vars=vars, color=color[1], label=label)
    if save:   
        path = tools.save_path(path)      
        plt.savefig(fname=os.path.join(path, prefix + '_' + win + '.png'))
    if show:
        plt.show()
    else:
        plt.close()






