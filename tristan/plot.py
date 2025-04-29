import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import dcmri as dc
import pydmr

from tristan import calc


sym = {
    'control': 'b-',
    'drug': 'r-',
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
    11: 'X',
    12: 'o',
    13: 'v',
    14: '^',
    15: '<',
    16: '>',
    17: 's',
    18: 'p',
    19: '*',
    20: 'x',
    21: 'd',
    22: 'X',
}


def color(index: float) -> str:
    cmap = plt.get_cmap("tab20")  
    rgba_color = cmap(index)  # Get RGBA tuple
    hex_color = "#{:02x}{:02x}{:02x}".format(
        int(rgba_color[0] * 255),
        int(rgba_color[1] * 255),
        int(rgba_color[2] * 255)
    )  # Convert to hex

    return hex_color



def _line_plot_ref(ax1, ax2, visits):

    output = pd.read_csv(
        os.path.join(os.getcwd(), 'tristan', 'reference.csv')
    )

    pivot = pd.pivot_table(output[output.visit=='control'], values='value', 
                           columns='parameter', index='subject')
    khe_ref = pivot.loc[:, 'khe']
    kbh_ref = pivot.loc[:, 'kbh']
    pivot = pd.pivot_table(output[output.visit=='drug'], values='value', 
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

    output = pd.read_csv(os.path.join(src, 'Analysis', 'parameters_rep.csv'))
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
    
    for i, s in enumerate(khe_ref.index):
        if s in khe_rif.index:
            x = [visits[0], visits[1]]
            khe = [khe_ref[s],khe_rif[s]]
            kbh = [kbh_ref[s],kbh_rif[s]]
        else:
            x = [visits[0]]
            khe = [khe_ref[s]]
            kbh = [kbh_ref[s]]    
        si = i/len(khe_ref.index.values)  
        ax1.plot(x, khe, '-', label=s, marker=mark[int(i+1)], 
                 markersize=markersize, color=color(si))
        ax2.plot(x, kbh, '-', label=s, marker=mark[int(i+1)], 
                 markersize=markersize, color=color(si))


def _effect_box_plots(src, ax):

    pars = ['khe', 'kbh']
    df = pd.read_csv(os.path.join(src, 'Analysis', 'effect_size.csv'))
    all_data = []
    for par in ['khe', 'kbh']:
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
    fig, (ax0, ax1, ax2) = plt.subplots(
        1, 3, width_ratios=[2, 4, 4], figsize=(8,3)
    )
    fig.subplots_adjust(wspace=0.5)
    _effect_box_plots(src, ax0)
    _line_plots(src, ax1, ax2, ylim=ylim, ref=ref)
    path = os.path.join(src, 'Figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(fname=os.path.join(path, '_effect_plot.png'))
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

    output_file = os.path.join(src, 'all_results')
    dmr = pydmr.read(output_file, format='table')
    df = pd.DataFrame(dmr['pars'], columns=['subject', 'visit', 'parameter', 'value'])

    df_ref = pd.read_csv(
        os.path.join(os.getcwd(), 'tristan', 'reference.csv')
    )

    # Add column and merge
    df['source'] = 'data'
    df_ref['source'] = 'reference'
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
      
    plt.savefig(fname=os.path.join(src, 'Figures', '_compare_to_ref.png'))
    plt.close()


def _ref_effect_box_plots(src, par, ax, ylim=[-100,0]):

    df = pd.read_csv(os.path.join(src, 'Analysis', 'effect_size.csv'))
    all_data = []
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
                        labels=['3-4hrs'])  # will be used to label x-ticks
    ax.set_ylim(ylim[0], ylim[1])
    ax.yaxis.grid(True)
    ax.tick_params(labelleft=False)
    if par == 'khe':
        ax.set_xticks([])

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    
    #ax.set_ylabel('Effect size (%)', fontsize=fontsize)


def _vart_effect_box_plots(src, par, ax, ylim=[-100,0]):

    df = pd.read_csv(os.path.join(src, 'Analysis', 'effect_size.csv'))
    all_data = []
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
    if par == 'khe':
        ax.set_xticks([])
    else:
        ax.set_xlabel('Total acquisition time (min)')

    # fill with colors
    for patch in bplot['boxes']:
        patch.set_facecolor('lightsteelblue')

    # adding horizontal grid line
    ax.yaxis.grid(True)
    ax.set_ylabel(par + ' effect size (%)', fontsize=fontsize)


def vart_effect_plot(src, src_2scan, ylim=None):
    if ylim is None:
        ylim = ([-100,-80], [-100,100])
    fig, ((ax0,ax1), (ax2,ax3)) = plt.subplots(
        2, 2, figsize=(8,4), width_ratios=[8,1],
    )
    fig.subplots_adjust(wspace=0.1)
    fig.suptitle('Effect sizes as a function of total acquisition time')
    _vart_effect_box_plots(src, 'khe', ax0, ylim=ylim[0])
    _vart_effect_box_plots(src, 'kbh', ax2, ylim=ylim[1])
    _ref_effect_box_plots(src_2scan, 'khe', ax1, ylim=ylim[0])
    _ref_effect_box_plots(src_2scan, 'kbh', ax3, ylim=ylim[1])
    path = os.path.join(src, 'Figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plt.savefig(fname=os.path.join(path, '_effect_plot.png'))
    plt.close()


def diurnal_k(src, ylim=[50,6]):

    output_file = os.path.join(src, 'all_results')
    dmr = pydmr.read(output_file, format='table')
    output = pd.DataFrame(dmr['pars'], columns=['subject', 'visit', 'parameter', 'value'])

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
        for i, s in enumerate(subjects):
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
                    si = i/len(subjects)
                    ax[visit+par].plot(
                        t, data_subj, '-', 
                        label=s, marker=mark[int(i+1)], 
                        markersize=markersize, color=color(si))
    path = os.path.join(src, 'Figures')
    if not os.path.exists(path):
        os.makedirs(path)
    plot_file = os.path.join(path, '_diurnal_function.png')
    plt.savefig(fname=plot_file)
    plt.close()


def create_bar_chart(resultsfolder, ylim={}):

    path = os.path.join(resultsfolder, 'Figures')
    if not os.path.exists(path):
        os.makedirs(path)

    output_file = os.path.join(resultsfolder, 'all_results')
    dmr = pydmr.read(output_file, format='table')
    output = pd.DataFrame(dmr['pars'], columns=['subject', 'visit', 'parameter', 'value'])

    output['group'] = calc.lookup(resultsfolder, output.parameter.values, 'group')
    output['description'] = calc.lookup(resultsfolder, output.parameter.values, 'description')
    output['unit'] = calc.lookup(resultsfolder, output.parameter.values, 'unit') 
    visits = output.visit[output.visit!='change (%)'].unique()

    # Create bar charts for each parameter
    subjects = output['subject'].unique()
    structures = output['group'].unique()
    for struct in structures:
        df_struct = output[output.group==struct]
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
                path, '_plot_' + par + '_' + struct + '.png'
            )
            plt.savefig(fname=plot_file)
            plt.close()