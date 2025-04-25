import os

from tristan import report, master

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    drug = 'ciclosporin'

    sourcepath = os.path.join(os.getcwd(), 'data', f'tristan_humans_healthy_{drug}.dmr')
    resultspath = os.path.join(outputpath, drug)

    master.run(
        sourcepath, 
        resultspath, 
        effect_range=([-100,200], [-100,500]),
        k_max=[100, 5],
        acq_times=[5,10,15,20],
        ref=True,
        compute=True,
    )
    report.build(
        resultspath, 
        drug,
        title = 'Two-compound study',
        subtitle = drug,
        subject = 'D2.10 - Internal report',
    )


if __name__ == '__main__':
    main()

