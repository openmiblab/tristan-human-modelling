import os
import shutil

import miblab

from methods import report, master


def main():

    drug = 'rifampicin'
    results = os.path.join(os.getcwd(), 'build', f'{drug}_clinical')
    datapath = os.path.join(os.getcwd(), 'data')

    data = miblab.zenodo_fetch(
        f'tristan_humans_patients_{drug}.dmr.zip', 
        datapath,
    )
    master.run(
        data, 
        results, 
        acq_times=[5,10,15,20],
        ref=True,
        compute=True,
    )
    report.all_results(
        results, 
        'report (complete)',
        title = 'Gothenburg patient study',
        subtitle = f'{drug} (all results)',
        subject = 'D2.07 - Internal report',
    )
    report.key_results(
        results, 
        'report (summary)',
        title = 'Gothenburg patient study',
        subtitle = f'{drug} (key results)',
        subject = 'D2.13 - Internal report',
    )
    shutil.rmtree(datapath)

    
if __name__ == '__main__':
    main()