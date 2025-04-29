import os

from tristan import onescan, twoscan, plot, tools, calc, tables


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
        path = os.path.join(resultspath, 'onescan')
        onescan.compute(dmr_file, path)

        # Twoscan
        path = os.path.join(resultspath, 'twoscan')
        twoscan.compute(dmr_file, path)

        # Variable time
        path = os.path.join(resultspath, 'onescan_vart')
        onescan.compute_vart(dmr_file, path, acq_times=acq_times)

    # Compute statistics
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        calc.derive_pars(path)
        calc.desc_stats(path)
        calc.pairwise_diff(path)
        calc.pairwise_ttest(path)

    # Create plots
    for exp in ['onescan','twoscan']:
        path = os.path.join(resultspath, exp)
        plot.create_bar_chart(path)
        plot.effect_plot(path, ylim=k_max, ref=ref)
    path = os.path.join(resultspath, 'twoscan')
    plot.diurnal_k(path, ylim=k_max)
    
    # Create tables
    for exp in ['onescan','twoscan']:
        tables.pairwise_diff(resultspath, exp)
        tables.pairwise_stats(resultspath, exp)
        tables.cases(resultspath, exp)
    
    # Variable acquisition time results
    path = os.path.join(resultspath, 'onescan_vart')
    calc.derive_vart_pars(path)
    plot.vart_effect_plot(
        path, os.path.join(resultspath, 'twoscan'), 
        ylim=effect_range,
    )

    if ref:
        # Compare to reference results
        path = os.path.join(resultspath, 'twoscan')
        calc.derive_pars(path, ref=True)
        calc.compare_to_ref(path)
        plot.compare_to_ref(path)
        tables.reference(resultspath, 'twoscan')
    else:
        # Generate reference data for future studies
        path = os.path.join(resultspath, 'twoscan')
        calc.derive_pars(path, ref=True)        