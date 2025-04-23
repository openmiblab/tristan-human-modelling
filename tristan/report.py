import os
import pandas as pd
import numpy as np
import pylatex as pl
from pylatex.utils import NoEscape

from tristan import calc


CONTROL = 'control'
DRUG = 'drug'

# # Needs harmonizing
# CONTROL = 'baseline'
# DRUG = 'rifampicin'

def section_diurnal(doc, results):
    with doc.create(pl.Section('Diurnal variation')):

        with doc.create(pl.Figure(position='h!')) as pic:
            im = os.path.join(results, 'twoscan', 'Plots', '_diurnal_function.png')
            pic.add_image(im, width='6in')
            pic.add_caption(
                "Intra-day changes in hepatocellular uptake (k_he, top row) "
                "and biliary excretion (k_bh, bottom row) of gadoxetate at "
                "for the control (left column) and treatment (right "
                "column). Full lines connect values taken in the same subject "
                "at the same day." )


def section_acqtime(doc, results):
    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Acquisition time')):

        with doc.create(pl.Figure(position='h!')) as pic:
            im = os.path.join(results, 'onescan_vart', 'Plots', '_effect_plot.png')
            pic.add_image(im, width='6in')
            pic.add_caption(
                "The effect of shortening the acquisition time on measured "
                "effect sizes. The figure shows the reference result with 2 "
                "scans, and then repeated results with shorter acquisition "
                "times ranging from 5min to 40min." )
            
            
def section_summary(doc, results, folder):

    with doc.create(pl.Section('Data summary')):

        with doc.create(pl.Figure(position='h!')) as pic:
            pic.append(pl.Command('centering'))
            im = os.path.join(results, folder, 'Plots', '_effect_plot.png')
            pic.add_image(im, width='7in')
            pic.add_caption(
                "Effect size (%) on hepatocellular uptake (k_he, left) and "
                "biliary excretion (k_bh, right) of gadoxetate. The boxplot "
                "shows median, interquartile range and 95 percent range. The "
                "line plots show individual values for hepatocellular uptake "
                "(k_he, middle) and biliary excretion (k_bh, right) of "
                "gadoxetate of the control (left of plot) and treatment "
                "(right of plot). Grey lines are healthy controls with "
                "rifampicin injection.")

        df = pd.read_csv(
            os.path.join(results, folder, 'Analysis', '_table_k_stats.csv')
        ).set_index('parameter')
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('l'+'c'*df.shape[1])) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption(
                "Effect size and absolute values of hepatocellular uptake "
                "(k_he) and biliary excretion (k_bh) of gadoxetate")
            

def section_reference(doc, results, folder):

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Comparison to reference values')):

        with doc.create(pl.Figure(position='h!')) as pic:
            im = os.path.join(results, folder, 'Plots', '_compare_to_ref.png')
            pic.add_image(im, width='6in')
            pic.add_caption(
                "Comparison to reference values in healthy volunteers "
                "treated with the drug and under control conditions.")

        doc.append(NoEscape('\\clearpage')) 
        with doc.create(pl.Subsection('Control')):
            
            df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats_ref.csv'))
            df = df[df.visit==CONTROL]
            df = df[df.group=='MRI - liver']
            df.drop(columns=['visit','group'], inplace=True)
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                    tab.add_hline()
                    tab.add_row([df.index.name] + list(df.columns))
                    tab.add_hline()
                    for row in df.index:
                        tab.add_row([row] + list(df.loc[row,:]))
                    tab.add_hline()
                table.add_caption(
                    "Results of a pairwise comparison testing for differences in "
                    "liver biomarkers between this study and healthy "
                    "reference data, under control conditions.")
                
            df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats_ref.csv'))
            df = df[df.visit==CONTROL]
            df = df[df.group=='MRI - aorta']
            df.drop(columns=['visit','group'], inplace=True)
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                    tab.add_hline()
                    tab.add_row([df.index.name] + list(df.columns))
                    tab.add_hline()
                    for row in df.index:
                        tab.add_row([row] + list(df.loc[row,:]))
                    tab.add_hline()
                table.add_caption(
                    "Results of a pairwise comparison testing for differences in "
                    "aorta biomarkers between this study and healthy "
                    "reference data, under control conditions.")

        doc.append(NoEscape('\\clearpage'))        
        with doc.create(pl.Subsection('Treatment')):
            
            df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats_ref.csv'))
            df = df[df.visit==DRUG]
            df = df[df.group=='MRI - liver']
            df.drop(columns=['visit','group'], inplace=True)
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                    tab.add_hline()
                    tab.add_row([df.index.name] + list(df.columns))
                    tab.add_hline()
                    for row in df.index:
                        tab.add_row([row] + list(df.loc[row,:]))
                    tab.add_hline()
                table.add_caption(
                    "Results of a pairwise comparison testing for differences in "
                    "liver biomarkers between this study and healthy "
                    "reference data, under treatment conditions.")
                
            df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats_ref.csv'))
            df = df[df.visit==DRUG]
            df = df[df.group=='MRI - aorta']
            df.drop(columns=['visit','group'], inplace=True)
            with doc.create(pl.Table(position='h!')) as table:
                table.append(pl.Command('centering'))
                with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                    tab.add_hline()
                    tab.add_row([df.index.name] + list(df.columns))
                    tab.add_hline()
                    for row in df.index:
                        tab.add_row([row] + list(df.loc[row,:]))
                    tab.add_hline()
                table.add_caption(
                    "Results of a pairwise comparison testing for differences in "
                    "aorta biomarkers between this study and healthy "
                    "reference data, under treatment conditions.")



def section_biomarkers(doc, results, folder):

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Liver biomarkers')):
        df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats2.csv')).set_index('Biomarker')
        df = df[df.group=='MRI - liver']
        df.drop(columns='group', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption(
                "Results of a pairwise comparison testing for differences in "
                "liver biomarkers between control and treatment. The "
                "results are ranked by their p-value, with most significant "
                "differences at the top of the list.")

        df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats1.csv')).set_index('Biomarker')
        df = df[df.group=='MRI - liver']
        df.drop(columns='group', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption(
                "Mean values along with their 95 percent confidence "
                "intervals for all liver biomarkers of the control and "
                "treatment visit. The last column shows the relative change "
                "at the treatment visit. The results are ranked by their "
                "p-value, with most significant differences at the top of the "
                "list.")

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Systemic biomarkers')):

        df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats2.csv')).set_index('Biomarker')
        df = df[df.group=='MRI - aorta']
        df.drop(columns='group', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption(
                "Results of a pairwise comparison testing for differences in "
                "systemic biomarkers between control and treatment visit. The "
                "results are ranked by their p-value, with most significant "
                "differences at the top of the list.")

        df = pd.read_csv(os.path.join(results, folder, 'Analysis', 'stats1.csv')).set_index('Biomarker')
        df = df[df.group=='MRI - aorta']
        df.drop(columns='group', inplace=True)
        with doc.create(pl.Table(position='h!')) as table:
            table.append(pl.Command('centering'))
            with table.create(pl.Tabular('ll'+'c'*(df.shape[1]-1))) as tab:
                tab.add_hline()
                tab.add_row([df.index.name] + list(df.columns))
                tab.add_hline()
                for row in df.index:
                    tab.add_row([row] + list(df.loc[row,:]))
                tab.add_hline()
            table.add_caption(
                "Mean values along with their 95 percent confidence "
                "intervals for all systemic biomarkers at the control and "
                "and treatment visit. The last column shows the relative "
                "change at the treatment visit. The results are ranked by "
                "their p-value, with most significant differences at the top "
                "of the list.")


def section_case_notes(doc, results, folder):

    doc.append(NoEscape('\\clearpage'))
    with doc.create(pl.Section('Case notes')):

        # Get data
        df = pd.read_csv(
            os.path.join(results, folder, 'Analysis', 'parameters_rep.csv')
        )
        src = os.path.join(results, folder)
        df['group'] = calc.lookup(src, df.parameter.values, 'group')
        df['name'] = calc.lookup(src, df.parameter.values, 'name')
        df['unit'] = calc.lookup(src, df.parameter.values, 'unit') 

        # Get effect sizes
        ef = pd.read_csv(
            os.path.join(results, folder, 'Analysis', 'effect_size.csv')
        )
        ef['group'] = calc.lookup(src, ef.parameter.values, 'group')
        ef['name'] = calc.lookup(src, ef.parameter.values, 'name')
        ef['visit'] = 'change (%)'
        ef['unit'] = calc.lookup(src, ef.parameter.values, 'unit')
        ef['stdev'] = np.nan

        # Concatenate
        ef = ef[['subject','visit','parameter','name','value','unit','stdev','group']]
        df = df[['subject','visit','parameter','name','value','unit','stdev','group']]
        df = pd.concat([df, ef])
        dfa = df[df.group=='MRI - aorta']
        df = df[df.group=='MRI - liver']

        for i, subject in enumerate(df.subject.unique()):
            if i>0:
                doc.append(NoEscape('\\clearpage'))
            subj = str(subject).zfill(3)
            with doc.create(pl.Subsection('Subject ' + subj)):

                dfs = df[df.subject==subject]
                dfas = dfa[dfa.subject==subject]
                dfr = dfs[dfs.visit==DRUG]

                # Images
                with doc.create(pl.Figure(position='h!')) as pic:
                    im = os.path.join(
                        results, folder, 'Figs', subj +'_'+CONTROL+'.png')
                    pic.add_image(im, width='4.5in')
                    pic.add_caption(
                        "Signal-time curves for subject "+subj+" at the "
                        "control visit.")

                if not dfr.empty:
                    with doc.create(pl.Figure(position='h!')) as pic:
                        im = os.path.join(
                            results, folder,  'Figs', subj +'_'+DRUG+'.png')
                        pic.add_image(im, width='4.5in')
                        pic.add_caption(
                            "Signal-time curves for subject "+subj+" at "
                            "the treatment visit.")

                # Tables
                doc.append(NoEscape('\\clearpage'))
                pivot = pd.pivot_table(dfs, values='value', columns='visit', 
                                       index=['name','unit'])
                cols = pivot.columns.tolist()
                if len(cols)>1:
                    pivot = pivot[[CONTROL,DRUG,'change (%)']]
                with doc.create(pl.Table(position='h!')) as table:
                    table.append(pl.Command('centering'))
                    tabular = pl.Tabular('ll'+'c'*pivot.shape[1])
                    with table.create(tabular) as tab:
                        rvals = ['Biomarker', 'Units'] + list(pivot.columns)
                        tab.add_hline()
                        tab.add_row(rvals)
                        tab.add_hline()
                        for row in pivot.index:
                            rvals = list(np.around(pivot.loc[row,:].values,2))
                            tab.add_row([row[0],row[1]] + rvals)
                        tab.add_hline()
                    table.add_caption("Values for liver of subject "+subj)

                pivot = pd.pivot_table(dfas, values='value', columns='visit', 
                                       index=['name','unit'])
                cols = pivot.columns.tolist()
                if len(cols)>1:
                    pivot = pivot[[CONTROL,DRUG,'change (%)']]
                with doc.create(pl.Table(position='h!')) as table:
                    table.append(pl.Command('centering'))
                    tabular = pl.Tabular('ll'+'c'*pivot.shape[1])
                    with table.create(tabular) as tab:
                        rvals = ['Biomarker', 'Units'] + list(pivot.columns)
                        tab.add_hline()
                        tab.add_row(rvals)
                        tab.add_hline()
                        for row in pivot.index:
                            rvals = list(np.around(pivot.loc[row,:].values,2))
                            tab.add_row([row[0],row[1]] + rvals)
                        tab.add_hline()
                    table.add_caption("Values for aorta of subject "+subj)


if __name__ == "__main__":
    generate()
