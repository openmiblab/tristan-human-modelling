from methods import twoscan, plot, calc, tables

def run(dmr_file, path, k_max=[100, 10]):

    # Compute all results
    twoscan.compute(dmr_file, path)

    # Compute statistics
    calc.effect_size(path)
    calc.descriptive_statistics(path)
    calc.averages(path)
    calc.pairwise_ttest(path)

    # Create plots
    plot.create_bar_chart(path)
    plot.effect_plot(path, ylim=k_max)
    plot.diurnal_k(path, ylim=k_max)
    
    # Create tables
    tables.averages(path)
    tables.pairwise_stats(path)
    tables.cases(path)
      