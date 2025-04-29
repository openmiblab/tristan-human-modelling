import os

import miblab

from tristan import onescan, twoscan

drug = 'controls'
results = os.path.join(os.getcwd(), 'results', drug)

def main(compute=True):

    data = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        os.path.join(os.getcwd(), 'data'),
    )

    if compute:

        # Onescan
        path = os.path.join(results, 'onescan')
        onescan.compute(data, path)

        # Twoscan
        path = os.path.join(results, 'twoscan')
        twoscan.compute(data, path)

        # Variable time
        path = os.path.join(results, 'onescan_vart')
        onescan.compute_vart(data, path)


if __name__ == '__main__':
    main()

