import os

import miblab

from tristan import report, master

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    drug = 'rifampicin'

    sourcepath = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )
    resultspath = os.path.join(outputpath, drug)

    master.run(
        sourcepath, 
        resultspath, 
        effect_range=([-100,200], [-100,500]),
        k_max=[100, 5],
        acq_times=[5,10,15,20,25,30,35,40],
        ref=False,
        compute=True,
    )
    report.build(
        resultspath, 
        drug,
        title = 'Experimental medicine study',
        subtitle = drug,
        subject = 'D2.13 - Internal report',
    )


if __name__ == '__main__':
    main()

