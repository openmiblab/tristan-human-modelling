import os

from tristan import onescan, twoscan 
import plot, tools, calc, report


root = os.path.abspath(os.sep)
datapath = os.path.join(root, 'Users', 'steve', 'Dropbox')
sourcepath = os.path.join(datapath, 'Data', 'tristan_gothenburg_final_data')
resultspath = os.path.join(datapath, 'Results', 'tristan_gothenburg_dev_final_data')

# Format data
onescan.format_data(sourcepath, os.path.join(resultspath, 'onescan'))
twoscan.format_data(sourcepath, os.path.join(resultspath, 'twoscan'))

# Calculate
tools.compute(onescan.fit_subj, os.path.join(resultspath, 'onescan'))
tools.compute(twoscan.fit_subj, os.path.join(resultspath, 'twoscan'))
tools.compute_vart(
    onescan.fit_subj, 
    os.path.join(resultspath, 'onescan_vart'),
    os.path.join(resultspath, 'onescan'),
    acq_times = [5,10,15,20])

# Summarise results
for exp in ['onescan','twoscan']:
    src = os.path.join(resultspath, exp)
    plot.create_bar_chart(src)
    calc.derive_pars(src)
    calc.ttest(src)
    calc.desc_stats(src)
    plot.effect_plot(src, ylim=[50,10], ref=True)

calc.derive_pars(os.path.join(resultspath, 'twoscan'), ref=True)
calc.compare_to_ref(resultspath)
plot.compare_to_ref(resultspath)

# Plot diurnal variations
plot.diurnal_k(os.path.join(resultspath, 'twoscan'), ylim=[50,10])

# # Variable acquisition time results
src = os.path.join(resultspath, 'onescan_vart')
plot.create_bar_chart(src)
calc.derive_vart_pars(src)
plot.vart_effect_plot(src, os.path.join(resultspath, 'twoscan'), 
                      ylim=([-100,200], [-100,500]))

# Create report
report.generate('gothenburg', resultspath)
report.generate('gothenburg', resultspath)