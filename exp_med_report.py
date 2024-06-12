import os
import pandas as pd
import numpy as np

import pylatex as pl
from pylatex.utils import NoEscape

#import miblab_report.report as report
from miblab_report import report

def generate(filename, results):

    print('Creating report..')

    # Cover and title pages
    report.setup(os.path.abspath(""), results)
    doc = pl.Document()
    doc.documentclass = pl.Command('documentclass',"epflreport")
    report.makecover(doc, 
            title = 'Experimental medicine study',
            subtitle = 'Results',
            subject = 'D2.13 - Internal report')
    report.titlepage(doc, results)
    doc.append(pl.NewPage())
    doc.append(NoEscape('\\tableofcontents'))
    doc.append(NoEscape('\\mainmatter'))

    # Two-scan results
    folder = 'twoscan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Two-scan results'))
    section_summary(doc, results, folder)
    section_biomarkers(doc, results, folder)
    section_case_notes(doc, results, folder)

    # One-scan results
    folder = 'onescan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'One-scan results'))
    section_summary(doc, results, folder)
    section_biomarkers(doc, results, folder)
    section_case_notes(doc, results, folder)

    # Secondary results
    #doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Secondary results'))
    with doc.create(pl.Section('Diurnal variation')):

        with doc.create(pl.Figure(position='h!')) as pic:
            im = os.path.join(results, 'twoscan', '_diurnal_function.png')
            pic.add_image(im, width='6in')
            pic.add_caption("Intra-day changes in hepatocellular uptake (k_he, top row) and biliary excretion (k_bh, bottom row) of Gadoxetate at baseline (left column) and after rifampicin (right column). Full lines connect values taken in the same subject at the same day." )

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Acquisition time')):

        with doc.create(pl.Figure(position='h!')) as pic:
            im = os.path.join(results, 'onescan_vart', '_effect_plot.png')
            pic.add_image(im, width='6in')
            pic.add_caption("The effect of shortening the acquisition time on measured effect sizes. The figure shows the reference result with 2 scans, and then repeated results with shorter acquisition times ranging from 5min to 40min." )

    report.create(doc, os.path.abspath(""), filename, results)


def section_summary(doc, results, folder):

    with doc.create(pl.Section('Data summary')):

        with doc.create(pl.Figure(position='h!')) as pic:
            pic.append(pl.Command('centering'))
            im = os.path.join(results, folder, '_effect_plot.png')
            pic.add_image(im, width='7in')
            pic.add_caption("Effect size (%) on hepatocellular uptake (k_he, left) and biliary excretion (k_bh, right) of Gadoxetate. The boxplot shows median, interquartile range and 95 percent range. The line plots show individual values for hepatocellular uptake (k_he, middle) and biliary excretion (k_bh, right) of Gadoxetate at baseline (left of plot) and after rifampicin (right of plot).")

        df = pd.read_pickle(os.path.join(results, folder, '_table_k_stats.pkl'))
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('l'+'c'*df.shape[1])) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption("Effect size and absolute values of hepatocellular uptake (k_he) and biliary excretion (k_bh) of Gadoxetate")

def section_biomarkers(doc, results, folder):

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Liver biomarkers')):

        df = pd.read_pickle(os.path.join(results, folder, 'stats2.pkl'))
        df = df[df.structure=='liver']
        df.drop(columns='structure', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption("Results of a pairwise comparison testing for differences in liver biomarkers between baseline visit and rifampicin. The results are ranked by their p-value, with most significant differences at the top of the list.")

        df = pd.read_pickle(os.path.join(results, folder, 'stats1.pkl'))
        df = df[df.structure=='liver']
        df.drop(columns='structure', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption("Mean values along with their 95 percent confidence intervals for all liver biomarkers at the baseline visit and after rifampicin. The last column shows the relative change induced by rifampicin. The results are ranked by their p-value, with most significant differences at the top of the list.")

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Systemic biomarkers')):

        df = pd.read_pickle(os.path.join(results, folder, 'stats2.pkl'))
        df = df[df.structure=='aorta']
        df.drop(columns='structure', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption("Results of a pairwise comparison testing for differences in systenic biomarkers between baseline visit and rifampicin. The results are ranked by their p-value, with most significant differences at the top of the list.")

        df = pd.read_pickle(os.path.join(results, folder, 'stats1.pkl'))
        df = df[df.structure=='aorta']
        df.drop(columns='structure', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption("Mean values along with their 95 percent confidence intervals for all systemic biomarkers at the baseline visit and after rifampicin. The last column shows the relative change induced by rifampicin. The results are ranked by their p-value, with most significant differences at the top of the list.")


def section_case_notes(doc, results, folder):

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Case notes')):
        df = pd.read_pickle(os.path.join(results, folder, 'parameters_rep.pkl'))
        dfa = df[df.structure=='aorta']
        df = df[df.structure=='liver']
        for i, subject in enumerate(df.subject.unique()):
            if i>0:
                doc.append(NoEscape('\\clearpage'))
            subj = str(subject).zfill(3)
            with doc.create(pl.Subsection('Subject ' + subj)):

                dfs = df[df.subject==subject]
                dfas = dfa[dfa.subject==subject]
                dfr = dfs[dfs.visit=='rifampicin']

                # Liver results
                with doc.create(pl.Figure(position='h!')) as pic:
                    im = os.path.join(results, folder, subj +'_baseline_Liver_all.png')
                    pic.add_image(im, width='5.5in')
                    pic.add_caption("Liver signal-time curves for subject "+subj+' at baseline.')
                if not dfr.empty:
                    with doc.create(pl.Figure(position='h!')) as pic:
                        im = os.path.join(results, folder,  subj +'_rifampicin_Liver_all.png')
                        pic.add_image(im, width='5.5in')
                        pic.add_caption("Liver signal-time curves for subject "+subj+' after rifampicin.')

                pivot = pd.pivot_table(dfs, values='value', columns='visit', index=['name','unit'])
                cols = pivot.columns.tolist()
                if len(cols)>1:
                    pivot = pivot[['baseline','rifampicin','change (%)']]
                with doc.create(pl.Table(position='h!')) as table:
                    table.append(pl.Command('centering'))
                    with table.create(pl.Tabular('ll'+'c'*pivot.shape[1])) as tab:
                        tab.add_hline()
                        tab.add_row(['Biomarker', 'Units'] + list(pivot.columns))
                        tab.add_hline()
                        for row in pivot.index:
                            tab.add_row([row[0],row[1]] + list(np.around(pivot.loc[row,:].values,2)))
                        tab.add_hline()
                    table.add_caption("Values for liver of subject "+subj)

                # Aorta results
                doc.append(NoEscape('\\clearpage'))
                with doc.create(pl.Figure(position='h!')) as pic:
                    im = os.path.join(results, folder,  subj +'_baseline_Aorta_all.png')
                    pic.add_image(im, width='5.5in')
                    pic.add_caption("Aorta signal-time curves for subject "+subj+' at baseline.')
                if not dfr.empty:
                    with doc.create(pl.Figure(position='h!')) as pic:
                        im = os.path.join(results, folder,  subj +'_rifampicin_Aorta_all.png')
                        pic.add_image(im, width='5.5in')
                        pic.add_caption("Aorta signal-time curves for subject "+subj+' after rifampicin.')

                pivot = pd.pivot_table(dfas, values='value', columns='visit', index=['name','unit'])
                cols = pivot.columns.tolist()
                if len(cols)>1:
                    pivot = pivot[['baseline','rifampicin','change (%)']]
                with doc.create(pl.Table(position='h!')) as table:
                    table.append(pl.Command('centering'))
                    with table.create(pl.Tabular('ll'+'c'*pivot.shape[1])) as tab:
                        tab.add_hline()
                        tab.add_row(['Biomarker', 'Units'] + list(pivot.columns))
                        tab.add_hline()
                        for row in pivot.index:
                            tab.add_row([row[0],row[1]] + list(np.around(pivot.loc[row,:].values,2)))
                        tab.add_hline()
                    table.add_caption("Values for aorta of subject "+subj)


if __name__ == "__main__":
    generate()
