import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


sym = {
    'baseline': 'b-',
    'rifampicin': 'r-',
}
mark = {
    1: 'o',
    2: 'v',
    3: '^',
    4: '<',
    5: '>',
    6: 's',
    7: 'p',
    8: '*',
    9: 'x',
    10: 'd',
}
clr = {
    1: 'tab:blue',
    2: 'tab:orange',
    3: 'tab:green',
    4: 'tab:red',
    5: 'tab:purple',
    6: 'tab:brown',
    7: 'tab:pink',
    8: 'tab:gray',
    9: 'tab:olive',
    10: 'tab:cyan',
}

def _line_plot_ref(ax1, ax2, visits):

    output = pd.read_pickle(
        os.path.join(os.getcwd(), 'tristan', 'reference.pkl'))

    pivot = pd.pivot_table(output[output.visit=='baseline'], values='value', 
                           columns='parameter', index='subject')
    khe_ref = pivot.loc[:, 'khe']
    kbh_ref = pivot.loc[:, 'kbh']
    pivot = pd.pivot_table(output[output.visit=='rifampicin'], values='value', 
                           columns='parameter', index='subject')
    khe_rif = pivot.loc[:, 'khe']
    kbh_rif = pivot.loc[:, 'kbh']

    for s in khe_ref.index:
        if s in khe_rif.index:
            x = [visits[0], visits[1]]
            khe = [khe_ref[s],khe_rif[s]]
            kbh = [kbh_ref[s],kbh_rif[s]]
        else:
            x = [visits[0]]
            khe = [khe_ref[s]]
            kbh = [kbh_ref[s]]      
        ax1.plot(x, khe, '-', label=s, marker='o', markersize=6, 
                 color='lightgrey')
        ax2.plot(x, kbh, '-', label=s, marker='o', markersize=6, 
                 color='lightgrey')
        

def _line_plots(src, ax1, ax2, ylim=[50,5], ref=False):

    output = pd.read_pickle(os.path.join(src, 'parameters_rep.pkl'))
    #output = pd.read_pickle(os.path.join(src, 'parameters_ext.pkl'))
    visits = output.visit[output.visit!='change (%)'].unique()

    fontsize=10
    markersize=6
    ax1.set_title('Hepatocellular uptake rate', fontsize=fontsize, pad=10)
    ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, ylim[0])
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=fontsize)
    ax2.set_title('Biliary excretion rate', fontsize=fontsize, pad=10)
    ax2.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax2.set_ylim(0, ylim[1])
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=fontsize)

    if ref:
        _line_plot_ref(ax1, ax2, visits)
    
    pivot = pd.pivot_table(output[output.visit==visits[0]], values='value', 
                           columns='parameter', index='subject')
    khe_ref = pivot.loc[:, 'khe']
    kbh_ref = pivot.loc[:, 'kbh']
    pivot = pd.pivot_table(output[output.visit==visits[1]], values='value', 
                           columns='parameter', index='subject')
    khe_rif = pivot.loc[:, 'khe']
    kbh_rif = pivot.loc[:, 'kbh']
    
    for s in khe_ref.index:
        if s in khe_rif.index:
            x = [visits[0], visits[1]]
            khe = [khe_ref[s],khe_rif[s]]
            kbh = [kbh_ref[s],kbh_rif[s]]
        else:
            x = [visits[0]]
            khe = [khe_ref[s]]
            kbh = [kbh_ref[s]]      
        ax1.plot(x, khe, '-', label=s, marker=mark[int(s)], 
                 markersize=markersize, color=clr[int(s)])
        ax2.plot(x, kbh, '-', label=s, marker=mark[int(s)], 
                 markersize=markersize, color=clr[int(s)])


def _effect_box_plots(src, ax):

    pars = ['khe', 'kbh']
    df = pd.read_pickle(os.path.join(src, 'parameters_rep.pkl'))
    #df = pd.read_pickle(os.path.join(src, 'parameters_ext.pkl'))
    all_data = []
    for par in ['khe effect size', 'kbh effect size']:
        data = df[df.parameter==par].value.values.tolist()
        all_data.append(data)

    linewidth = 1.0
    fontsize=10
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    # box plot
    boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
    medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
    whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
    capprops = dict(linestyle='-', linewidth=linewidth, color='black')
    flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                  markeredgecolor='black', markeredgewidth=linewidth)
    bplot = ax.boxplot(all_data,
                        whis = [2.5,97.5],
                        capprops=capprops,
                        flierprops=flierprops,
                        whiskerprops=whiskerprops,
                        medianprops=medianprops,
                        boxprops=boxprops,
                        widths=0.3,
                        vert=True,  # vertical box alignment
                        patch_artist=True,  # fill with color
                        labels=pars)  # will be used to label x-ticks
    ax.set_xticklabels(labels=pars, fontsize=fontsize)
    ax.set_yticklabels(labels=ax.get_yticklabels(), fontsize=fontsize)

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    ax.yaxis.grid(True)
    ax.set_ylabel('Effect size (%)', fontsize=fontsize)


def effect_plot(src, ylim=[50,5], ref=False):
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, width_ratios=[2, 4, 4], 
                                        figsize=(8,3))
    fig.subplots_adjust(wspace=0.5)
    _effect_box_plots(src, ax0)
    _line_plots(src, ax1, ax2, ylim=ylim, ref=ref)
    plt.savefig(fname=os.path.join(src, '_effect_plot.png'))
    plt.close()


def _compare_to_ref_box_plot(ax, all_data, ylabel=None, title=None, ylim=None):

    pars = ['reference', 'study']
    
    linewidth = 1.0
    fontsize=10
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    # box plot
    boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
    medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
    whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
    capprops = dict(linestyle='-', linewidth=linewidth, color='black')
    flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                  markeredgecolor='black', markeredgewidth=linewidth)
    ax.set_ylim(ylim)
    bplot = ax.boxplot(all_data,
                        whis = [2.5,97.5],
                        capprops=capprops,
                        flierprops=flierprops,
                        whiskerprops=whiskerprops,
                        medianprops=medianprops,
                        boxprops=boxprops,
                        widths=0.3,
                        vert=True,  # vertical box alignment
                        patch_artist=True,  # fill with color
                        labels=pars)  # will be used to label x-ticks
    ax.set_xticklabels(labels=pars, fontsize=fontsize)
    ax.set_yticklabels(labels=ax.get_yticklabels(), fontsize=fontsize)

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    ax.yaxis.grid(True)
    
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=fontsize)
    if title is not None:
        ax.set_title(title, fontsize=12)
    


def compare_to_ref(src):

    df = pd.read_pickle(os.path.join(src, 'twoscan', 'reference.pkl'))
    df_ref = pd.read_pickle(
        os.path.join(os.getcwd(), 'tristan', 'reference.pkl'))
    
    # This becomes obsolete afte renaming exp_med visits
    df_ref.replace('baseline', 'control', inplace=True)
    df_ref.replace('rifampicin', 'drug', inplace=True)

    # Add column and merge
    df_ref['source'] = 'reference'
    df['source'] = 'data'
    df = pd.concat((df_ref, df))

    # Plot boxes
    fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2, 2, figsize=(5,5), 
                                               width_ratios=[1,1])
    #fig.subplots_adjust(wspace=0.1)
    #fig.suptitle('Comparison against reference data')

    ax = {
        'khe': {
            'control': ax0,
            'drug': ax1,
        },
        'kbh': {
            'control': ax2,
            'drug': ax3,
        },
    }
    ylim = {
        'khe': [0, 60],
        'kbh': [0, 6],
    }

    for par in ['khe', 'kbh']:
        dfp = df[df.parameter==par]
        for visit in ['control', 'drug']:
            dfv = dfp[dfp.visit==visit]
            x = dfv[dfv.source=='reference'].value.values.tolist()
            y = dfv[dfv.source=='data'].value.values.tolist()
            ylabel = par + ' (mL/min/100mL)' if visit=='control' else None
            title = visit if par=='khe' else None
            _compare_to_ref_box_plot(ax[par][visit], [x,y], ylabel, title, ylim[par])
      
    plt.savefig(fname=os.path.join(src, 'twoscan', '_compare_to_ref.png'))
    plt.close()


def _ref_effect_box_plots(src, biom, ax, ylim=[-100,0]):

    df = pd.read_pickle(os.path.join(src, 'parameters_rep.pkl'))
    all_data = []
    data = df[df.parameter==biom + ' effect size'].value.values.tolist()
    all_data.append(data)

    linewidth = 1.0
    fontsize=10
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    # box plot
    boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
    medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
    whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
    capprops = dict(linestyle='-', linewidth=linewidth, color='black')
    flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                  markeredgecolor='black', markeredgewidth=linewidth)
    bplot = ax.boxplot(all_data,
                        whis = [2.5,97.5],
                        capprops=capprops,
                        flierprops=flierprops,
                        whiskerprops=whiskerprops,
                        medianprops=medianprops,
                        boxprops=boxprops,
                        widths=0.3,
                        vert=True,  # vertical box alignment
                        patch_artist=True,  # fill with color
                        labels=['3-4hrs'])  # will be used to label x-ticks
    ax.set_ylim(ylim[0], ylim[1])
    ax.yaxis.grid(True)
    ax.tick_params(labelleft=False)
    if biom == 'khe':
        ax.set_xticks([])

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    
    #ax.set_ylabel('Effect size (%)', fontsize=fontsize)


def _vart_effect_box_plots(src, biom, ax, ylim=[-100,0]):

    df = pd.read_pickle(os.path.join(src, 'parameters_rep.pkl'))
    all_data = []
    par = biom + ' effect size'
    labels = df.tacq.unique()
    for tacq in labels:
        data = df[(df.parameter==par) & (df.tacq==tacq)]
        data = data.value.values.tolist()
        all_data.append(data)

    linewidth = 1.0
    fontsize=10
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    # box plot
    boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
    medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
    whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
    capprops = dict(linestyle='-', linewidth=linewidth, color='black')
    flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                  markeredgecolor='black', markeredgewidth=linewidth)
    bplot = ax.boxplot(all_data,
                        whis = [2.5,97.5],
                        capprops=capprops,
                        flierprops=flierprops,
                        whiskerprops=whiskerprops,
                        medianprops=medianprops,
                        boxprops=boxprops,
                        widths=0.3,
                        vert=True,  # vertical box alignment
                        patch_artist=True,  # fill with color
                        labels=labels)  # will be used to label x-ticks
    ax.set_ylim(ylim[0], ylim[1])
    if biom == 'khe':
        ax.set_xticks([])
    else:
        ax.set_xlabel('Total acquisition time (min)')

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    ax.yaxis.grid(True)
    ax.set_ylabel(biom + ' effect (%)', fontsize=fontsize)


def vart_effect_plot(src, src_2scan, ylim=None):
    if ylim is None:
        ylim = ([-100,-80], [-100,100])
    fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(2, 2, figsize=(8,4), 
                                               width_ratios=[8,1])
    fig.subplots_adjust(wspace=0.1)
    fig.suptitle('Effect sizes as a function of total acquisition time')
    _vart_effect_box_plots(src, 'khe', ax0, ylim=ylim[0])
    _vart_effect_box_plots(src, 'kbh', ax2, ylim=ylim[1])
    _ref_effect_box_plots(src_2scan, 'khe', ax1, ylim=ylim[0])
    _ref_effect_box_plots(src_2scan, 'kbh', ax3, ylim=ylim[1])
    plt.savefig(fname=os.path.join(src, '_effect_plot.png'))
    plt.close()


def _max_line_plots(output_file, ax1, ax2):

    resultsfolder = os.path.dirname(output_file)
    output = pd.read_pickle(os.path.join(resultsfolder, 'parameters_ext.pkl'))
    visits = output.visit[output.visit!='change (%)'].unique()

    fontsize=12
    markersize=6
    ax1.set_title('Hepatocellular uptake rate', fontsize=fontsize, pad=10)
    ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, 40)
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=fontsize)
    ax2.set_title('Biliary excretion rate', fontsize=fontsize, pad=10)
    ax2.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax2.set_ylim(0, 3.5)
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=fontsize)

    pivot = pd.pivot_table(output[output.visit==visits[0]], values='value', 
                           columns='parameter', index='subject')
    khe_ref = pivot.loc[:, 'khe_max']
    kbh_ref = pivot.loc[:, 'kbh_max']
    pivot = pd.pivot_table(output[output.visit==visits[1]], values='value', 
                           columns='parameter', index='subject')
    khe_rif = pivot.loc[:, 'khe_min']
    kbh_rif = pivot.loc[:, 'kbh_min']
    
    for s in khe_ref.index:
        if s in khe_rif.index:
            x = visits
            khe = [khe_ref[s],khe_rif[s]]
            kbh = [kbh_ref[s],kbh_rif[s]]
        else:
            x = [visits[0]]
            khe = [khe_ref[s]]
            kbh = [kbh_ref[s]]            
        ax1.plot(x, khe, '-', label=s, marker=mark[int(s)], 
                 markersize=markersize, color=clr[int(s)])
        ax2.plot(x, kbh, '-', label=s, marker=mark[int(s)], 
                 markersize=markersize, color=clr[int(s)])
    #ax1.legend(loc='upper center', ncol=5, prop={'size': 14})
    #ax2.legend(loc='upper center', ncol=5, prop={'size': 14})


def _max_effect_box_plots(output_file, ax):

    resultsfolder = os.path.dirname(output_file)
    output = pd.read_pickle(os.path.join(resultsfolder,
                                          'parameters_ext.pkl'))
    khe = output[output.parameter=='khe max effect size 1'].value.values.tolist()
    kbh = output[output.parameter=='kbh max effect size 1'].value.values.tolist()
    all_data = [khe, kbh]
    parameters = ['khe', 'kbh']
    
    # Create box plot
    linewidth = 1.0
    fontsize=10
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(linewidth)

    # box plot
    boxprops = dict(linestyle='-', linewidth=linewidth, color='black')
    medianprops = dict(linestyle='-', linewidth=linewidth, color='black')
    whiskerprops = dict(linestyle='-', linewidth=linewidth, color='black')
    capprops = dict(linestyle='-', linewidth=linewidth, color='black')
    flierprops = dict(marker='o', markerfacecolor='white', markersize=6,
                  markeredgecolor='black', markeredgewidth=linewidth)
    bplot = ax.boxplot(all_data,
                        whis = [2.5,97.5],
                        capprops=capprops,
                        flierprops=flierprops,
                        whiskerprops=whiskerprops,
                        medianprops=medianprops,
                        boxprops=boxprops,
                        widths=0.3,
                        vert=True,  # vertical box alignment
                        patch_artist=True,  # fill with color
                        labels=parameters)  # will be used to label x-ticks
    
    ax.set_xticklabels(labels=parameters, fontsize=fontsize)
    ax.set_yticklabels(labels=ax.get_yticklabels(), fontsize=fontsize)

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('white')

    # adding horizontal grid line
    ax.yaxis.grid(True)
    ax.set_ylabel('Maximum effect size (%)', fontsize=fontsize)
    #ax.set_ylim(-100, 20)


def max_effect_plot(output_file):
    fig, (ax0, ax1, ax2) = plt.subplots(1, 3, width_ratios=[2, 4, 4], 
                                        figsize=(8,3))
    fig.subplots_adjust(wspace=0.5)
    _max_effect_box_plots(output_file, ax0)
    _max_line_plots(output_file, ax1, ax2)
    resultsfolder = os.path.dirname(output_file)
    plt.savefig(fname=os.path.join(resultsfolder, '_max_effect_plot.png'))
    plt.close()


def tmax_effect(output_file, pdf):
    # Read parameter file
    output = pd.read_pkl(output_file)
    output = output[output.structure=='liver']
    # Set up figure
    fontsize=12
    titlesize=14
    markersize=2
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(8.27,11.69))
    fig.subplots_adjust(
                    left=0.1,
                    right=0.95,
                    bottom=0.2,
                    top = 0.80, 
                    wspace=0.3,
                    #hspace=1,
                    )
    title =( 
        "Effect of changing acquisition time \non hepatocellular uptake "
        "(khe, top row) and biliary excretion (k_bh, bottom row) of "
        "Gadoxetate \n at baseline (left column) and after administration "
        "of rifampicin (right column). \n Full lines connect values taken in "
        "the same subject." )
    fig.suptitle(title, fontsize=12)
    ax = {
        'baselinekhe': ax1,
        'rifampicinkhe': ax2,
        'baselinekbh': ax3,
        'rifampicinkbh': ax4,
    }
    ax1.set_title('Baseline', fontsize=titlesize)
    ax1.set_xlabel('Scan time (mins)', fontsize=fontsize)
    ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, 50)
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=fontsize)
    ax2.set_title('Rifampicin', fontsize=titlesize)
    ax2.set_xlabel('Scan time (mins)', fontsize=fontsize)
    ax2.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax2.set_ylim(0, 3)
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=fontsize)
    #ax3.set_title('Baseline', fontsize=titlesize)
    ax3.set_xlabel('Scan time (mins)', fontsize=fontsize)
    ax3.set_ylabel('k_bh (mL/min/100mL)', fontsize=fontsize)
    ax3.set_ylim(0, 4)
    ax3.tick_params(axis='x', labelsize=fontsize)
    ax3.tick_params(axis='y', labelsize=fontsize)
    #ax4.set_title('Rifampicin', fontsize=titlesize)
    ax4.set_xlabel('Scan time (mins)', fontsize=fontsize)
    ax4.set_ylabel('k_bh (mL/min/100mL)', fontsize=fontsize)
    ax4.set_ylim(0, 4)
    ax4.tick_params(axis='x', labelsize=fontsize)
    ax4.tick_params(axis='y', labelsize=fontsize)
    
    # Create plots
    for visit in output.visit.unique():
        for subject in output.subject.unique():
            df = output[output.visit==visit]
            df = df[df.subject==subject]
            if not df.empty:
                df = pd.pivot(df, values='value', columns='parameter', 
                              index='tmax') 
                for par in ['khe', 'kbh']:
                    ax[visit+par].plot(
                        df.index.values/60, df[par].values, 'k-', 
                        label=subject, marker=mark[int(subject)], 
                        markersize=markersize)

    # Export results
    resultsfolder = os.path.dirname(output_file)
    plot_file = os.path.join(resultsfolder, '_tmax.png')
    plt.savefig(fname=plot_file)
    pdf.savefig()
    #plt.show()
    plt.close()
                

def diurnal_k(src, ylim=[50,6]):
    # output = pd.read_pickle(os.path.join(src, 'parameters_ext.pkl'))
    # visits = output.visit[output.visit!='change (%)'].unique()

    output = pd.read_pickle(os.path.join(src, 'parameters.pkl'))
    visits = output.visit.unique()

    fontsize=10
    titlesize=12
    markersize=6
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2,2,figsize=(8,8))
    fig.subplots_adjust(
                    left=0.1,
                    right=0.9,
                    bottom=0.1,
                    top = 0.9, 
                    wspace=0.3,
                    #hspace=1,
                    )
    ax = {
        visits[0]+'khe': ax1,
        visits[1]+'khe': ax2,
        visits[0]+'kbh': ax3,
        visits[1]+'kbh': ax4,
    }
    ax1.set_title(visits[0], fontsize=titlesize)
    ax1.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax1.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, ylim[0])
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=fontsize)
    ax2.set_title(visits[1], fontsize=titlesize)
    ax2.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax2.set_ylabel('khe (mL/min/100mL)', fontsize=fontsize)
    ax2.set_ylim(0, ylim[0])
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=fontsize)
    #ax3.set_title('Baseline', fontsize=titlesize)
    ax3.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax3.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax3.set_ylim(0, ylim[1])
    ax3.tick_params(axis='x', labelsize=fontsize)
    ax3.tick_params(axis='y', labelsize=fontsize)
    #ax4.set_title('Rifampicin', fontsize=titlesize)
    ax4.set_xlabel('Time of day (hrs)', fontsize=fontsize)
    ax4.set_ylabel('kbh (mL/min/100mL)', fontsize=fontsize)
    ax4.set_ylim(0, ylim[1])
    ax4.tick_params(axis='x', labelsize=fontsize)
    ax4.tick_params(axis='y', labelsize=fontsize)

    # Create box plots
    for visit in visits:
        df_visit = output[output.visit==visit]
        subjects = df_visit.subject.unique()
        for s in subjects:
            df_subj = df_visit[df_visit.subject==s]
            for par in ['khe', 'kbh']:
                data_subj = []
                for p in [par+'_i', par+'_f']:
                    df_par = df_subj[df_subj.parameter==p]
                    if not df_par.empty:
                        v = df_par.value.values[0]
                        data_subj.append(v)
                t = []
                for p in ['t0', 't3']:
                    df_par = df_subj[df_subj.parameter==p]
                    if not df_par.empty:
                        v = df_par.value.values[0]
                        t.append(v)
                if len(data_subj) == 2:
                    ax[visit+par].plot(
                        t, data_subj, '-', 
                        label=s, marker=mark[int(s)], 
                        markersize=markersize, color=clr[int(s)])

    plot_file = os.path.join(src, '_diurnal_function.png')
    plt.savefig(fname=plot_file)
    plt.close()



def line_plot_extracellular(output_file):
    resultsfolder = os.path.dirname(output_file)
    output = pd.read_pickle(output_file)
    fontsize=20
    titlesize=30
    fig, (ax1, ax2) = plt.subplots(1,2,figsize=(16,8))
    fig.tight_layout(pad=10.0)
    fig.suptitle('Extracellular drug effects', fontsize=titlesize)
    ax1.set_title('Volume fraction', fontsize=24)
    ax1.set_ylabel('ve (mL/100mL)', fontsize=fontsize)
    ax1.set_ylim(0, 60)
    ax1.tick_params(axis='x', labelsize=fontsize)
    ax1.tick_params(axis='y', labelsize=18)
    ax2.set_title('Transit time', fontsize=24)
    ax2.set_ylabel('Te (min)', fontsize=fontsize)
    ax2.set_ylim(0, 1)
    ax2.tick_params(axis='x', labelsize=fontsize)
    ax2.tick_params(axis='y', labelsize=18)

    ax = {
        've': ax1,
        'Te': ax2,
    }
    subjects = output['subject'].unique()
    visits = output.visit[output.visit!='change (%)'].unique()
    structures = output['structure'].unique()
    for struct in structures:
        df_struct = output[output.structure==struct]
        for par in ['ve', 'Te']:
            df_par = df_struct[df_struct.parameter==par]
            for s in subjects:
                df_subj = df_par[df_par.subject==s]
                x = []
                y = []
                for visit in visits:
                    df_visit = df_subj[df_subj.visit==visit]
                    if not df_visit.empty:
                        v = df_visit['value'].values[0]
                        x.append(visit)
                        y.append(v)
                ax[par].plot(x, y, 'k-', label=s, marker=mark[int(s)], 
                             markersize=12)
    #ax1.legend(loc='upper center', ncol=5, prop={'size': 14})
    #ax2.legend(loc='upper center', ncol=5, prop={'size': 14})
    plot_file = os.path.join(resultsfolder, '_lineplot_extracellular.png')
    #plt.show()
    plt.savefig(fname=plot_file)
    plt.close()


def create_box_plot(output_file, ylim={}):

    resultsfolder = os.path.dirname(output_file)
    output = pd.read_pickle(output_file)

    # Create box plots for each parameter
    subjects = output['subject'].unique()
    visits = output.visit[output.visit!='change (%)'].unique()
    structures = output['structure'].unique()
    for struct in structures:
        df_struct = output[output.structure==struct]
        for par in df_struct['parameter'].unique():
            df = df_struct[df_struct.parameter==par]
            all_data = []
            for visit in visits:
                df_visit = df[df.visit==visit]
                data_visit = []
                for s in subjects:
                    df_visit_subj = df_visit[df_visit.subject==s]
                    if not df_visit_subj.empty:
                        val = df_visit_subj['value'].values[0]
                    data_visit.append(val)
                all_data.append(data_visit)

            fig, ax = plt.subplots(layout='constrained')

            # notch shape box plot
            bplot = ax.boxplot(all_data,
                                #notch=True,  # notch shape
                                vert=True,  # vertical box alignment
                                patch_artist=True,  # fill with color
                                labels=visits)  # will be used to label x-ticks
            ax.set_title('Drug effect on ' + par)

            # fill with colors
            colors = ['slateblue', 'coral']
            for patch, color in zip(bplot['boxes'], colors):
                patch.set_facecolor(color)

            # adding horizontal grid line
            ax.yaxis.grid(True)
            ax.set_xlabel('Visit')
            units = df.unit.unique()[0]
            ax.set_ylabel(par + ' (' + str(units) + ')')
            #ax.legend(loc='upper left', ncols=2)
            if par in ylim:
                ax.set_ylim(ylim[par][0], ylim[par][1])

            plot_file = os.path.join(
                resultsfolder, '_boxplot_' + struct + '_' + par + '.png')
            plt.savefig(fname=plot_file)
            plt.close()


def create_bar_chart(resultsfolder, ylim={}):

    output_file = os.path.join(resultsfolder, 'parameters.pkl')
    output = pd.read_pickle(output_file)
    visits = output.visit[output.visit!='change (%)'].unique()

    # Create bar charts for each parameter
    subjects = output['subject'].unique()
    structures = output['structure'].unique()
    for struct in structures:
        df_struct = output[output.structure==struct]
        for par in df_struct['parameter'].unique():
            if par == 'Kbh':
                # For some reason the file with kbh is not written
                # when Kbh is already saved. Same for Khe - not written 
                # because khe is written first. Putting in this hack as 
                # Kbh is not of interest.
                continue
            df = df_struct[df_struct.parameter==par]
            bar_chart = {}
            for visit in visits:
                df_visit = df[df.visit==visit]
                bar_chart[visit] = []
                for s in subjects:
                    df_visit_subj = df_visit[df_visit.subject==s]
                    if df_visit_subj.empty:
                        val = np.nan
                    else:
                        val = df_visit_subj['value'].values[0]
                    bar_chart[visit].append(val)
            x = np.arange(len(subjects))  # the label locations
            width = 0.25  # the width of the bars
            multiplier = 0

            fig, ax = plt.subplots(layout='constrained')
            colors = {visits[0]:'slateblue', visits[1]:'coral'}
            for attribute, measurement in bar_chart.items():
                offset = width * multiplier
                if measurement != np.nan:
                    rects = ax.bar(
                        x + offset, measurement, width, label=attribute, 
                        color=colors[attribute])
                    ax.bar_label(rects, padding=3)
                multiplier += 1

            # Add some text for labels, title and custom x-axis tick labels, etc.
            units = df.unit.unique()[0]
            ax.set_ylabel(par + ' (' + str(units) + ')')
            ax.set_title('Values per visit ' + par)
            ax.set_xticks(x + width, subjects)
            ax.legend(loc='upper left', ncols=2)
            if par in ylim:
                ax.set_ylim(ylim[par][0], ylim[par][1])
            
            plot_file = os.path.join(
                resultsfolder, '_plot_' + par + '_' + struct + '.png')
            plt.savefig(fname=plot_file)
            plt.close()



# def drug_effect(output_file, pars=None, name='', ylim=[-100,100]):

#     resultsfolder = os.path.dirname(output_file)
#     output = pd.read_pickle(output_file)

#     # Create box plots for aorta and liver
#     subjects = output['subject'].unique()
#     structures = output['structure'].unique()
#     for struct in structures:
#         df_struct = output[output.structure==struct]
#         all_data = []
#         if pars is None:
#             parameters = df_struct['parameter'].unique()
#         else:
#             parameters = pars
#         for par in parameters:
#             df_par = df_struct[df_struct.parameter==par]
#             data_par = []
#             for s in subjects:
#                 df_subj = df_par[df_par.subject==s]
#                 df_baseline_subj = df_subj[df_subj.visit=='control']
#                 df_rifampicin_subj = df_subj[df_subj.visit=='drug'] 
#                 if not df_baseline_subj.empty and not df_rifampicin_subj.empty:
#                     v0 = df_baseline_subj['value'].values[0]
#                     v1 = df_rifampicin_subj['value'].values[0]
#                     if v0 != 0:
#                         data_par.append(100*(v1-v0)/v0)
#             all_data.append(data_par)

#         fig, ax = plt.subplots(layout='constrained')

#         # notch shape box plot
#         bplot = ax.boxplot(all_data,
#                             vert=True,  # vertical box alignment
#                             patch_artist=True,  # fill with color
#                             labels=parameters)  # will be used to label x-ticks
#         ax.set_title('Drug effect on ' + struct)

#         # fill with colors
#         for patch in bplot['boxes']:
#             patch.set_facecolor('blue')

#         # adding horizontal grid line
#         ax.yaxis.grid(True)
#         ax.set_xlabel('Parameter')
#         ax.set_ylabel('Rifampicin effect (%)')
#         ax.set_ylim(ylim[0], ylim[1])

#         plot_file = os.path.join(resultsfolder, 
#                                  '_drug_effect_' +name+ '_' + struct + '.png')
#         plt.savefig(fname=plot_file)
#         plt.close()



if __name__ == "__main__":
    tmax_effect()