"""
Replicate all results from TRISTAN studies in humans.
"""

import tristan_controls
import tristan_rifampicin
import tristan_rifampicin_clinical
import tristan_metformin
import tristan_ciclosporin

tristan_controls.main()
tristan_rifampicin_clinical.main()
tristan_rifampicin.main()
tristan_metformin.main()
tristan_ciclosporin.main()
