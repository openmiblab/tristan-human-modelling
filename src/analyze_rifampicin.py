import os
import shutil

import miblab

from methods import report, master_primary


def main():

    drug = 'rifampicin'
    dataset = f'tristan_humans_healthy_{drug}.dmr.zip'

    # Define folders
    results = os.path.join(os.getcwd(), 'build', drug)
    datapath = os.path.join(os.getcwd(), 'data')
    
    # Get the data
    data = miblab.zenodo_fetch(dataset, datapath)

    # Compute results
    master_primary.run(data, results)

    # Report results
    report.primary_results(
        results, 
        'report (analysis)',
        title = 'Leeds pilot study',
        subtitle = f'{drug} (analysis)',
        subject = 'Internal report',
    )

    # Cleanup temporary folder
    shutil.rmtree(datapath)


if __name__ == '__main__':
    main()

