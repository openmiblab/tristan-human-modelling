import os

from tristan import onescan, twoscan, plot, tools, calc
import exp_med_report

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    sourcepath = os.path.join(os.getcwd(), 'data', 'rifampicin.dmr')
    resultspath = os.path.join(outputpath, 'tristan_exp_med')

    # Onescan
    path = os.path.join(resultspath, 'onescan')
    onescan.compute(sourcepath, path)

    # Twoscan
    path = os.path.join(resultspath, 'twoscan')
    twoscan.compute(sourcepath, path)

    # Variable time
    path = os.path.join(resultspath, 'onescan_vart')
    onescan.compute_vart(sourcepath, path)

    # Compute statistics
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        tools.build_master(path)
        calc.derive_pars(path)
        calc.desc_stats(path)
        calc.ttest(path)

    # Create plots
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        plot.create_bar_chart(path)
        plot.effect_plot(path, ylim=[50,5])
    plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,5])

    # Variable acquisition time results
    path = os.path.join(resultspath, 'onescan_vart')
    tools.build_master(path, vart=True)
    calc.derive_vart_pars(path)
    plot.vart_effect_plot(path, os.path.join(resultspath, 'twoscan'))

    # Generate reference data for future studies
    calc.derive_pars(os.path.join(resultspath, 'twoscan'), ref=True)

    # Create report
    exp_med_report.build(resultspath)
    exp_med_report.build(resultspath)


if __name__ == '__main__':
    main()

