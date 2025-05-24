import os
import shutil

from studies import (
    controls,
    rifampicin,
    rifampicin_clinical,
    metformin,
    ciclosporin,
)

compute=True

# Compute all available studies
controls.main(compute)
rifampicin_clinical.main(compute)
rifampicin.main(compute)
metformin.main(compute)
ciclosporin.main(compute)

# Clean up temporary folder
data = os.path.join(os.getcwd(), 'data')
shutil.rmtree(data)