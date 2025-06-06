"""
Run the analysis of the TRISTAN project on a new drug.
"""

import os

from methods import report, master_primary


def main():

    drug = 'newdrug'
    data = 'path\to\newdrug.dmr'

    # Define folders
    results = os.path.join(os.getcwd(), 'build', drug)
    
    # Compute results
    master_primary.run(data, results)

    # Report results
    report.primary_results(
        results, 
        'report (analysis)',
        title = f'{drug} pilot study',
        subtitle = f'{drug} (analysis)',
        subject='Internal report',
    )


if __name__ == '__main__':
    main()

