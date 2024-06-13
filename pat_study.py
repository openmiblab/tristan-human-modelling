import os

from tristan import onescan, twoscan 
import plot, calc, pat_study_report

root = os.path.abspath(os.sep)
datapath = os.path.join(root, 'Users', 'steve', 'Dropbox')
sourcepath = os.path.join(datapath, 'Data', 'tristan_gothenburg')
resultspath = os.path.join(datapath, 'Results', 'tristan_gothenburg_dev')

# Calculate
# onescan.main(sourcepath, resultspath)
# twoscan.main(sourcepath, resultspath)

# Summarise results
for exp in ['onescan','twoscan']:
    src = os.path.join(resultspath, exp)
    calc.derive_pars(src)
    calc.ttest(src)
    calc.desc_stats(src)
    plot.effect_plot(src, ylim=[50,10])

# Plot diurnal variations
plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,10])

# Create report
pat_study_report.generate('report_gothenburg', resultspath)
pat_study_report.generate('report_gothenburg', resultspath)