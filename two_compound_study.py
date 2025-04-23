import os

from tristan import onescan, twoscan, plot, tools, calc
import two_compound_study_report

def main():

    root = os.path.abspath(os.sep)
    outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

    for drug in ['ciclosporin', 'metformin']:

        sourcepath = os.path.join(os.getcwd(), 'data', drug + '.dmr')
        resultspath = os.path.join(outputpath, 'tristan_two_compounds', drug)

        # Onescan
        path = os.path.join(resultspath, 'onescan')
        onescan.compute(sourcepath, path)

        # Twoscan
        path = os.path.join(resultspath, 'twoscan')
        twoscan.compute(sourcepath, path)

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
            plot.effect_plot(path, ylim=[100,10], ref=True)
        plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[100,10])
        plot.compare_to_ref(resultspath)

        # Create report
        two_compound_study_report.build(drug, resultspath)
        two_compound_study_report.build(drug, resultspath)


if __name__ == '__main__':
    main()

