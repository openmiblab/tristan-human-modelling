import os

from methods import onescan, twoscan, plot, calc, tables

ONESCAN = 'results (one scan)'
TWOSCAN = 'results (two scans)'
VART = 'results (one scan - variable tacq)'


def run(
        dmr_file, 
        resultspath, 
        effect_range=([-100,200], [-100,500]),
        k_max = [100, 10],
        acq_times=[5,10,15,20,25,30,35,40],
        ref=True,
        compute=True,
    ):

    if compute:

        # Onescan
        path = os.path.join(resultspath, ONESCAN)
        onescan.compute(dmr_file, path)

        # Twoscan
        path = os.path.join(resultspath, TWOSCAN)
        twoscan.compute(dmr_file, path)

        # Variable time
        path = os.path.join(resultspath, VART)
        onescan.compute_vart(dmr_file, path, acq_times=acq_times)

    # Compute statistics
    for exp in [ONESCAN, TWOSCAN]:
        path = os.path.join(resultspath, exp)
        calc.effect_size(path)
        calc.descriptive_statistics(path)
        calc.averages(path)
        calc.pairwise_ttest(path)

    # Create plots
    for exp in [ONESCAN, TWOSCAN]:
        path = os.path.join(resultspath, exp)
        plot.create_bar_chart(path)
        plot.effect_plot(path, ylim=k_max, ref=ref)
    path = os.path.join(resultspath, TWOSCAN)
    plot.diurnal_k(path, ylim=k_max)
    
    # Create tables
    for exp in [ONESCAN, TWOSCAN]:
        path = os.path.join(resultspath, exp)
        tables.averages(path)
        tables.pairwise_stats(path)
        tables.cases(path)
    
    # Variable acquisition time results
    path = os.path.join(resultspath, VART)
    calc.derive_vart_pars(path)
    plot.vart_effect_plot(
        path, os.path.join(resultspath, TWOSCAN), 
        ylim=effect_range,
    )

    path = os.path.join(resultspath, TWOSCAN)
    if ref:
        # Compare to reference results
        calc.effect_size(path, ref=True)
        calc.compare_to_ref(path)
        plot.compare_to_ref(path)
        tables.reference(path)
    else:
        # Generate reference data for future studies
        calc.effect_size(path, ref=True)        