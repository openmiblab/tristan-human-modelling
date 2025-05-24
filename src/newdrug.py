from studies import unknown


# The full path to the dmr file with your data
#data='path\to\newdrug.dmr',
data = 'C:\\Users\\md1spsx\\Downloads\\newdrug.dmr'

# Name of the output folder in build
drug = 'unknown_drug'

# Title of the pdf report
title = 'Unknown drug study'

# Subtitle of the pdf report
subtitle = 'Final results'

# Subject matter of the pdf report
subject = 'Internal report'


unknown.main(data, drug, title, subtitle, subject)

