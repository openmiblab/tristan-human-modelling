import os
import time

# Note - small bug in old version causes some smaller differences with previous results
# The issue was in the handling of the 3d MOLLI in 2-scan aorta model.
# It was transformed with S02 but predicted with S01. 

#Move parametrization back from Th to Kbh?

from tristan import aorta as fit_aorta_1scan
from tristan import liver as fit_liver_1scan
from tristan import aorta_2scan as fit_aorta_2scan
from tristan import liver_2scan as fit_liver_2scan

import calc
import plot
import pat_study_report

tstart = time.time()

root = os.path.abspath(os.sep)
datapath = os.path.join(root, 'Users', 'steve', 'Dropbox')
sourcepath = os.path.join(datapath, 'Data', 'tristan_gothenburg')
resultspath = os.path.join(datapath, 'Results', 'tristan_gothenburg_dev')

# # Calculate
# fit_aorta_1scan.main(sourcepath, resultspath)
fit_liver_1scan.main(sourcepath, resultspath)
# fit_aorta_2scan.main(sourcepath, resultspath)
fit_liver_2scan.main(sourcepath, resultspath)

# Aorta 2-scan results
src = os.path.join(resultspath, 'aorta_2scan')
calc.derive_aorta_pars(src)
calc.ttest(src, 'parameters_ext.pkl')

# Liver 2-scan results
src = os.path.join(resultspath, 'liver_2scan')
calc.derive_liver_pars(src)
calc.report_pars(src)
calc.desc_stats(src)
plot.effect_plot(src, ylim=[50,7])
plot.diurnal_k(src, ylim=[50,7])
# calc.ttest(src, 'parameters_rep.pkl')

# Aorta 1-scan results
src = os.path.join(resultspath, 'aorta_1scan')
calc.derive_aorta_pars(src)
#calc.ttest(src, 'parameters_ext.pkl')

# Liver 1-scan results
src = os.path.join(resultspath, 'liver_1scan')
calc.derive_liver_pars(src) 
calc.desc_stats(src)
plot.effect_plot(src, ylim=[50,7])
#calc.ttest(src, 'parameters_ext.pkl')

print('Computation time (mins)', (time.time()-tstart)/60)

pat_study_report.generate('report_gothernburg', resultspath)