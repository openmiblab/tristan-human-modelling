import os

from tristan import onescan, twoscan 
import plot, tools, calc, report

root = os.path.abspath(os.sep)
outputpath = os.path.join(root, 'Users', 'md1spsx', 'Documents', 'Results')

#for drug in ['ciclosporin', 'metformin']:
for drug in ['metformin']:

    sourcepath = os.path.join(os.getcwd(), 'data', 'tristan_two_compounds', drug)
    resultspath = os.path.join(outputpath, 'tristan_two_compounds', drug)

    # Format data
    onescan.format_data(sourcepath, os.path.join(resultspath, 'onescan'))
    twoscan.format_data(sourcepath, os.path.join(resultspath, 'twoscan'))

    # Calculate
    tools.compute(onescan.fit_subj, os.path.join(resultspath, 'onescan'))
    tools.compute(twoscan.fit_subj, os.path.join(resultspath, 'twoscan'))

    # Summarise results
    for exp in ['onescan','twoscan']:
        src = os.path.join(resultspath, exp)
        plot.create_bar_chart(src)
        calc.derive_pars(src)
        calc.ttest(src)
        calc.desc_stats(src)
        plot.effect_plot(src, ylim=[100,10], ref=True)

    # Compare to reference data
    calc.derive_pars(os.path.join(resultspath, 'twoscan'), ref=True)
    calc.compare_to_ref(resultspath)
    plot.compare_to_ref(resultspath)

    # Plot diurnal variations
    plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[100,10])

    # Create report
    report.generate(drug, resultspath)
    report.generate(drug, resultspath)

