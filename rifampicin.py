import os

import miblab

from tristan import report, master

drug = 'rifampicin'
results = os.path.join(os.getcwd(), 'results', drug)

def main(compute=True):

    data = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )
    master.run(
        data, 
        results, 
        ref=False,
        compute=compute,
    )
    report.all_results(
        results, 
        drug + '_all_results',
        title = 'Experimental medicine study',
        subtitle = drug + ' (all results)',
        subject = 'D2.13 - Internal report',
    )
    report.key_results(
        results, 
        drug + '_key_results',
        title = 'Experimental medicine study',
        subtitle = drug + ' (key results)',
        subject = 'D2.13 - Internal report',
    )


if __name__ == '__main__':
    main()

