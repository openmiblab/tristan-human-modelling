import os
import shutil

import miblab

from methods import onescan, twoscan


def main():

    drug = 'controls'
    results = os.path.join(os.getcwd(), 'build', drug)
    datapath = os.path.join(os.getcwd(), 'data')

    data = miblab.zenodo_fetch(
        f'tristan_humans_healthy_{drug}.dmr.zip', 
        datapath,
    )

    # One scan
    path = os.path.join(results, 'results (one scan)')
    onescan.compute(data, path)

    # Two scan
    path = os.path.join(results, 'results (two scans)')
    twoscan.compute(data, path)

    # Variable time
    path = os.path.join(results, 'results (one scan - variable tacq)')
    onescan.compute_vart(data, path)

    shutil.rmtree(datapath)


if __name__ == '__main__':
    main()

