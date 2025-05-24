import os

from methods import report, master


def main(
        data='path\to\dmrfile', 
        drug='unknown_drug', 
        title='Unknown drug study',
        subtitle='Final results',
        subject='Internal report',
        compute=True,
    ):

    results = os.path.join(os.getcwd(), 'build', drug)

    master.run(
        data, 
        results, 
        ref=True,
        compute=compute,
    )
    report.all_results(
        results, 
        'report (complete)',
        title=title,
        subtitle=f'{subtitle} (all results)',
        subject=subject,
    )
    report.key_results(
        results, 
        'report (summary)',
        title=title,
        subtitle=f'{subtitle} (key results)',
        subject=subject,
    )


if __name__ == '__main__':
    main()

