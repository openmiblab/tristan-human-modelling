import os

from tristan import onescan, twoscan
import plot, tools, calc, report

root = os.path.abspath(os.sep)
datapath = os.path.join(root, 'Users', 'steve', 'Dropbox')
sourcepath = os.path.join(datapath, 'Data', 'tristan_exp_med')
resultspath = os.path.join(datapath, 'Results', 'tristan_exp_med_dev')

# Format data
# onescan.format_data(sourcepath, os.path.join(resultspath, 'onescan'))
twoscan.format_data(sourcepath, os.path.join(resultspath, 'twoscan'))

# Calculate
# tools.compute(onescan.fit_subj, os.path.join(resultspath, 'onescan'))
tools.compute(twoscan.fit_subj, os.path.join(resultspath, 'twoscan'))
# tools.compute_vart(
#     onescan.fit_subj, 
#     os.path.join(resultspath, 'onescan_vart'),
#     os.path.join(resultspath, 'onescan'))

# Summarise results
#for exp in ['onescan','twoscan']:
for exp in ['twoscan']:
    src = os.path.join(resultspath, exp)
    plot.create_bar_chart(src)
    calc.derive_pars(src)
    calc.ttest(src)
    calc.desc_stats(src)
    plot.effect_plot(src, ylim=[50,5])

# Plot diurnal variations
# plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,5])

# # Variable acquisition time results
# src = os.path.join(resultspath, 'onescan_vart')
# plot.create_bar_chart(src)
# calc.derive_vart_pars(src)
# plot.vart_effect_plot(src, os.path.join(resultspath, 'twoscan'))

# # Create report
# report.generate('exp_med', resultspath)
# report.generate('exp_med', resultspath)

# Generate reference data for future studies
calc.derive_pars(os.path.join(resultspath, 'twoscan'), ref=True)