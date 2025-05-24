import os

import miblab

from methods import onescan, twoscan

ONESCAN = 'results (one scan)'
TWOSCAN = 'results (two scans)'
VART = 'results (one scan - variable tacq)'

drug = 'controls'
results = os.path.join(os.getcwd(), 'build', drug)

def main(compute=True):

    data = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )

    if compute:

        # Onescan
        path = os.path.join(results, ONESCAN)
        onescan.compute(data, path)

        # Twoscan
        path = os.path.join(results, TWOSCAN)
        twoscan.compute(data, path)

        # Variable time
        path = os.path.join(results, VART)
        onescan.compute_vart(data, path)


if __name__ == '__main__':
    main()

