import os

import miblab

from tristan import report, master

drug = 'rifampicin'
results = os.path.join(os.getcwd(), 'results', f'{drug}_clinical')

def main(compute=True):

    data = miblab.zenodo_fetch(
        f'tristan_humans_patients_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )
    master.run(
        data, 
        results, 
        acq_times=[5,10,15,20],
        ref=True,
        compute=compute,
    )
    report.all_results(
        results, 
        f'{drug}_clinical' + '_all_results',
        title = 'Gothenburg patient study',
        subtitle = drug + ' (all results)',
        subject = 'D2.07 - Internal report',
    )
    report.key_results(
        results, 
        f'{drug}_clinical' + '_key_results',
        title = 'Gothenburg patient study',
        subtitle = drug + ' (key results)',
        subject = 'D2.13 - Internal report',
    )

    
if __name__ == '__main__':
    main()