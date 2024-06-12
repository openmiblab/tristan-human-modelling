import os
import pandas as pd
import numpy as np

import pylatex as pl
from pylatex.utils import NoEscape

import miblab_report.report as report


def generate(filename, results, subtitle='DCE-MRI results'):

    print('Creating report..')

    report.setup(os.path.abspath(""), results)

    doc = pl.Document()
    doc.documentclass = pl.Command('documentclass',"epflreport")

    report.makecover(doc, 
            title = 'Ghent dog kidney study',
            subtitle = subtitle,
            subject = 'Internal report',
            author = 'Steven Sourbron')
    report.titlepage(doc, results)
    # TOC page
    doc.append(pl.NewPage())
    doc.append(NoEscape('\\tableofcontents'))
    doc.append(NoEscape('\\mainmatter'))

    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Kidney biomarkers'))

    df = pd.read_pickle(os.path.join(results, 'kidneys', 'stats2.pkl'))
    df = df.reset_index()
    with doc.create(pl.Table(position='h!')) as table:
        table.append(pl.Command('centering'))
        with table.create(pl.Tabular('lll'+'c'*(df.shape[1]-3))) as tab:
            tab.add_hline()
            tab.add_row(list(df.columns))
            tab.add_hline()
            for row in df.index:
                tab.add_row(list(df.loc[row,:]))
            tab.add_hline()
        table.add_caption("Results of a pairwise comparison testing for differences in kidney biomarkers between baseline visit and followup. The results are ranked by their p-value, with most significant differences at the top of the list.")

    df = pd.read_pickle(os.path.join(results, 'kidneys', 'stats1.pkl'))
    df = df.reset_index()
    with doc.create(pl.Table(position='h!')) as table:
        table.append(pl.Command('centering'))
        with table.create(pl.Tabular('lll'+'c'*(df.shape[1]-3))) as tab:
            tab.add_hline()
            tab.add_row(list(df.columns))
            tab.add_hline()
            for row in df.index:
                tab.add_row(list(df.loc[row,:]))
            tab.add_hline()
        table.add_caption("Mean values along with their 95 percent confidence intervals for all kidney biomarkers at the baseline visit and at follow-up. The results are ranked by their p-value, with most significant differences at the top of the list.")

    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Aorta biomarkers'))

    df = pd.read_pickle(os.path.join(results, 'aorta', 'stats2.pkl'))
    df = df.reset_index()
    with doc.create(pl.Table(position='h!')) as table:
        table.append(pl.Command('centering'))
        with table.create(pl.Tabular('lll'+'c'*(df.shape[1]-3))) as tab:
            tab.add_hline()
            tab.add_row(list(df.columns))
            tab.add_hline()
            for row in df.index:
                tab.add_row(list(df.loc[row,:]))
            tab.add_hline()
        table.add_caption("Results of a pair wise comparison testing for differences in systemic biomarkers between baseline visit and follow-up. The results are ranked by their p-value, with most significant differences at the top of the list.")

    df = pd.read_pickle(os.path.join(results, 'aorta', 'stats1.pkl'))
    df = df.reset_index()
    with doc.create(pl.Table(position='h!')) as table:
        table.append(pl.Command('centering'))
        with table.create(pl.Tabular('lll'+'c'*(df.shape[1]-3))) as tab:
            tab.add_hline()
            tab.add_row(list(df.columns))
            tab.add_hline()
            for row in df.index:
                tab.add_row(list(df.loc[row,:]))
            tab.add_hline()
        table.add_caption("Mean values along with their 95 percent confidence intervals for all systemic biomarkers at the baseline visit and after rifampicin. The results are ranked by their p-value, with most significant differences at the top of the list.")


    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Case notes'))

    df = pd.read_pickle(os.path.join(results, 'kidneys', 'parameters.pkl'))
    dfa = pd.read_pickle(os.path.join(results, 'aorta', 'parameters_ext.pkl'))
    dfa = dfa.astype({'subject':'int32'})

    for i, subject in enumerate(df.subject.unique()):
        if i>0:
            doc.append(NoEscape('\\clearpage'))
        subj = str(subject)
        with doc.create(pl.Subsection('Subject ' + subj)):
            
            dfs = df[df.subject==subject]
            dfas = dfa[dfa.subject==int(subject)]
            
            # # Left Kidney plots
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'kidneys',  'Dog' + subj + '.1 LK_kidney_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Left kidney signal-time curves for subject "+subj+' at baseline.')
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'kidneys',  'Dog' + subj + '.2 LK_kidney_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Left kidney signal-time curves for subject "+subj+' at visit 2.')

            # Right Kidney plots
            doc.append(NoEscape('\\clearpage'))
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'kidneys',  'Dog' + subj + '.1 RK_kidney_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Right kidney signal-time curves for subject "+subj+' at baseline.')
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'kidneys',  'Dog' + subj + '.2 RK_kidney_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Right kidney signal-time curves for subject "+subj+' at visit 2.')

            # Numerical results for both kidneys
            doc.append(NoEscape('\\clearpage'))
            pivot = pd.pivot_table(dfs, values='value', columns='visit', index=['name','structure','unit'])
            cols = pivot.columns.tolist()
            if len(cols)>1:
                pivot = pivot[['1','2']]
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('lll'+'c'*pivot.shape[1])) as tab:
                    tab.add_hline()
                    tab.add_row(['Biomarker', 'Kidney', 'Units'] + list(pivot.columns))
                    tab.add_hline()
                    for row in pivot.index:
                        tab.add_row([row[0],row[1],row[2]] + list(np.around(pivot.loc[row,:].values,2)))
                    tab.add_hline()
                table.add_caption("Values for left kidney of subject "+subj)

            # Aorta plots
            doc.append(NoEscape('\\clearpage'))
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'aorta', 'Dog' + subj + '.1 AIF_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Aorta signal-time curves for subject "+subj+' at baseline.')
            with doc.create(pl.Figure(position='h!')) as pic:
                im = os.path.join(results, 'aorta',  'Dog' + subj + '.2 AIF_all.png')
                pic.add_image(im, width='5.5in')
                pic.add_caption("Aorta signal-time curves for subject "+subj+' at visit2.')

            # Numerical results for aorta
            pivot = pd.pivot_table(dfas, values='value', columns='visit', index=['name','structure','unit'])
            cols = pivot.columns.tolist()
            if len(cols)>1:
                pivot = pivot[['1','2']]
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('lll'+'c'*pivot.shape[1])) as tab:
                    tab.add_hline()
                    tab.add_row(['Biomarker', 'ROI', 'Units'] + list(pivot.columns))
                    tab.add_hline()
                    for row in pivot.index:
                        tab.add_row([row[0],row[1],row[2]] + list(np.around(pivot.loc[row,:].values,2)))
                    tab.add_hline()
                table.add_caption("Values for aorta of subject "+subj)

    report.create(doc, os.path.abspath(""), filename, results)


if __name__ == "__main__":
    generate()
