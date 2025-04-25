import os

from tristan import report, master

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    drug = 'rifampicin'

    sourcepath = os.path.join(os.getcwd(), 'data', f'tristan_humans_patients_{drug}.dmr')
    resultspath = os.path.join(outputpath, f'{drug}_clinical')

    master.run(
        sourcepath, 
        resultspath, 
        effect_range=([-100,200], [-100,500]),
        k_max=[50, 5],
        acq_times=[5,10,15,20],
        ref=True,
    )
    report.build(
        resultspath, 
        f'{drug}_clinical',
        title = 'Gothenburg patient study',
        subtitle = drug,
        subject = 'D2.07 - Internal report',
    )

    
if __name__ == '__main__':
    main()