import os

from tristan import onescan, twoscan, plot, tools, calc
import pat_study_report

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    sourcepath = os.path.join(os.getcwd(), 'data', 'rifampicin-patients.dmr')
    resultspath = os.path.join(outputpath, 'tristan_gothenburg')

    # Onescan
    path = os.path.join(resultspath, 'onescan')
    onescan.compute(sourcepath, path)

    # Twoscan
    path = os.path.join(resultspath, 'twoscan')
    twoscan.compute(sourcepath, path)

    # Variable time
    path = os.path.join(resultspath, 'onescan_vart')
    onescan.compute_vart(sourcepath, path, acq_times = [5,10,15,20])

    # Compute statistics
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        tools.build_master(path)
        calc.derive_pars(path)
        calc.desc_stats(path)
        calc.ttest(path)
    path = os.path.join(resultspath, 'twoscan')
    calc.derive_pars(path, ref=True)
    calc.compare_to_ref(path)

    # Create plots
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        plot.create_bar_chart(path)
        plot.effect_plot(path, ylim=[50,10], ref=True)
    plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,10])
    plot.compare_to_ref(resultspath)

    # Variable acquisition time results
    path = os.path.join(resultspath, 'onescan_vart')
    tools.build_master(path, vart=True)
    calc.derive_vart_pars(path)
    plot.vart_effect_plot(
        path, os.path.join(resultspath, 'twoscan'), 
        ylim=([-100,200], [-100,500])
    )

    # Create report
    pat_study_report.build(resultspath)
    pat_study_report.build(resultspath)

    
if __name__ == '__main__':
    main()