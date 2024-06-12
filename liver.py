import os
import time

import numpy as np
import pandas as pd
import dcmri as dc

import plot
from tristan import data
from dcmods import tools, fig

class Liver(dc.Liver):
    def pars(self, aorta:dc.Aorta):

        aorta.tmax = aorta.BAT+180*60
        self.t, self.ca = aorta.conc()
        self.ca = self.ca/(1-self.Hct)
        t, R1 = self.relax()
        t, C = self.conc()
        AUC_DR1 = np.trapz(R1-self.R10, self.t)
        AUC_C = np.trapz(C, self.t) 

        aorta.tmax = aorta.BAT+35*60
        self.t, self.ca = aorta.conc()
        self.ca = self.ca/(1-self.Hct)
        t, R1 = self.relax()
        t, C = self.conc()
        AUC35_DR1 = np.trapz(R1-self.R10, self.t)
        AUC35_C = np.trapz(C, self.t)  

        pars = super().pars()
        pars['AUC_R1']=['AUC for DR1 (0-inf)', AUC_DR1, '',0]
        pars['AUC_C']=['AUC for C (0-inf)', 1000*AUC_C, 'mM*sec',0]
        pars['AUC35_R1']=['AUC for DR1 (0-35min)', AUC35_DR1, '',0]
        pars['AUC35_C']=['AUC for C (0-35min)', 1000*AUC35_C, 'mM*sec',0]
        return pars   

def fit(data, path, name, aortapars):

    # Format data
    xdata = np.array(data['time1'])
    xvalid = np.array(data['liver_valid1'])
    ydata = np.array(data['liver1'])
    xcheck = np.array([0, data['T1time2']])
    R1 = [1000.0/data['T1liver1'], 
          1000.0/data['T1liver2']]
    xdata = xdata[xvalid==1]
    ydata = ydata[xvalid==1]

    # Get aorta model
    dt = 0.5
    aorta = dc.Aorta(
        dt = dt,
        tmax = max(xdata)+xdata[1]+dt,
        weight = data['weight'],
        agent = 'gadoxetate',
        dose = data['dose1'],
        rate = 1,
        field_strength = 3.0,
        TR = 3.71/1000.0,
        FA = 15,
        R10 = 1000.0/data['T1aorta1'],
        t0 = data['baseline'],
    )

    # Read aorta parameters
    df = tools.read_csv(aortapars)
    for par in aorta.free:
        setattr(aorta, par, df.at[par,"value"])

    # Generate input concentrations
    t, cb = aorta.conc()

    # For liver model
    liver = Liver(cb=cb,
        dt = dt,
        field_strength = 3.0,
        TR = 3.71/1000.0,
        FA = 15,
        agent = 'gadoxetate',
        Hct = 0.45,
        R10 = 1000.0/data['T1liver1'],
        t0 = data['baseline'],
        vol = data['liver_volume'],
    )

    # Fit model to data
    loss0 = liver.cost(xdata, ydata)
    print('Goodness of fit (initial): ', loss0)
    liver.train(xdata, ydata, xtol=1e-4, ftol=1e-4, verbose=2)
    loss1 = liver.cost(xdata, ydata)
    print('Goodness of fit (improvement, %): ', 100*(loss0-loss1)/loss0)

    # Export data
    # ycheck = [
    #     dc.signal_ss(R1[0], liver.S0, liver.TR, liver.FA, R10=R1[0]),
    #     dc.signal_ss(R1[1], liver.S0, liver.TR, liver.FA, R10=R1[0])]
    ycheck = [
        dc.signal_ss(R1[0], liver.S0, liver.TR, liver.FA),
        dc.signal_ss(R1[1], liver.S0, liver.TR, liver.FA)]
    fig.liver1scan_all(aorta.BAT, 
            liver.predict, liver.conc, xdata, ydata, 
            xcheck=xcheck, ycheck=ycheck, path=path, prefix=name+'_Liver')

    pars = liver.pars(aorta)
    tools.to_csv(liver, os.path.join(path, name + '.csv'), pars)
    pars = tools.to_tristan_units(pars)
    return tools.to_df(pars)


def main(datapath, results):

    start = time.time()
    resultspath = os.path.join(results, 'liver_1scan')
    aortaresults = os.path.join(results, 'aorta_1scan') 

    output = None
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            print('Fitting liver of ', visit, subj)
            subj_data = data.read(subj)
            name = s[:3] + '_' + visit
            aortapars = os.path.join(aortaresults, name + '.csv')
            liver_pars = fit(subj_data, resultspath, name, aortapars)
            liver_pars['subject'] = s[:3]
            liver_pars['visit'] = visit
            liver_pars['structure'] = 'liver'
            if output is None:
                output = liver_pars
            else:
                output = pd.concat([output, liver_pars])

    # Format output and save
    output = output.reindex(columns=['subject','visit','structure','name','value','unit','stdev'])
    output['parameter'] = output.index
    output.to_csv(os.path.join(resultspath, 'parameters.csv'), index=False)
    output.to_pickle(os.path.join(resultspath, 'parameters.pkl'))
    plot.create_bar_chart(os.path.join(resultspath, 'parameters.pkl'))
    
    print('Calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()