import os
import time

import pandas as pd

import plot, tools
from tristan import data, onescan


def main(datapath, results):

    start = time.time()
    resultspath = os.path.join(results, 'onescan_vart') 
    resultspath = tools.save_path(resultspath)

    output = None
    for visit in [f.name for f in os.scandir(datapath) if f.is_dir()]:
        visitdatapath = os.path.join(datapath, visit)
        for s in os.listdir(visitdatapath):
            subj = os.path.join(visitdatapath, s)
            subj_data = data.read(subj)
            for tacq in [5,10,15,20,25,30,35,40]:
                print('Fitting data of ', visit, subj, tacq)
                name = s[:3] + '_' + visit + '_' + str(tacq)
                pars = onescan.fit(subj_data, resultspath, name, tacq*60)
                pars['subject'] = s[:3]
                pars['visit'] = visit
                structure = []
                for p in pars.index.values:
                    if p in ['BAT','CO','Thl','Dhl','To',
                            'Eb','Eo','Teb','BAT2','Tc',
                            'AUC_R1b','AUC_Cb',
                            'AUC35_R1b','AUC35_Cb']:
                        structure.append('aorta')
                    else:
                        structure.append('liver')
                pars['structure'] = structure
                pars['tacq'] = tacq
                if output is None:
                    output = pars
                else:
                    output = pd.concat([output, pars])

    # Format output and save
    output = output.reindex(columns=['subject','visit','structure','name','value','unit','stdev','tacq'])
    output['parameter'] = output.index
    output.to_csv(os.path.join(resultspath, 'parameters.csv'), index=False)
    output.to_pickle(os.path.join(resultspath, 'parameters.pkl'))
    plot.create_bar_chart(os.path.join(resultspath, 'parameters.pkl'))
    
    print('Calculation time (mins): ', (time.time()-start)/60)


if __name__ == "__main__":
    main()