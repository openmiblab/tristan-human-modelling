import os

import miblab

from methods import report, master

drug = 'ciclosporin'
results = os.path.join(os.getcwd(), 'build', drug)

def main(compute=True):

    data = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )
    master.run(
        data, 
        results, 
        ref=True,
        compute=compute,
    )
    report.all_results(
        results, 
        'report (complete)',
        title = 'Sheffield two-compound study',
        subtitle = drug + ' (all results)',
        subject = 'D2.10 - Internal report',
    )
    report.key_results(
        results, 
        'report (summary)',
        title = 'Sheffield two-compound study',
        subtitle = drug + ' (key results)',
        subject = 'D2.13 - Internal report',
    )


if __name__ == '__main__':
    main()

