import os

from tristan import onescan, twoscan, onescan_vart 
import plot, calc, exp_med_report

root = os.path.abspath(os.sep)
datapath = os.path.join(root, 'Users', 'steve', 'Dropbox')
sourcepath = os.path.join(datapath, 'Data', 'tristan_exp_med')
resultspath = os.path.join(datapath, 'Results', 'tristan_exp_med_dev')

# Calculate
onescan.main(sourcepath, resultspath)
twoscan.main(sourcepath, resultspath)
onescan_vart.main(sourcepath, resultspath)

# Summarise results
for exp in ['onescan','twoscan']:
    src = os.path.join(resultspath, exp)
    calc.derive_pars(src)
    calc.ttest(src)
    calc.desc_stats(src)
    plot.effect_plot(src, ylim=[50,5])

# Plot diurnal variations
plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,5])

# Plot variable acquisition time results
src = os.path.join(resultspath, 'onescan_vart')
calc.derive_vart_pars(src)
plot.vart_effect_plot(src, os.path.join(resultspath, 'twoscan'))

# Creat report
exp_med_report.generate('report_exp_med', resultspath)
exp_med_report.generate('report_exp_med', resultspath)