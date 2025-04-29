import os
import pandas as pd
import miblab


def key_results(
        resultspath, 
        filename,
        title = 'Liver-mediated DDI study',
        subtitle = 'Key results',
        subject = 'Internal report',
    ):

    print('Creating report..')

    # Cover and title pages
    doc = miblab.Report(
        resultspath,
        filename,
        title = title,
        subtitle = subtitle,
        subject = subject,
    )

    folder = 'twoscan'
    doc.chapter('Key results')
    section_summary(doc, resultspath, folder)
    section_biomarkers(doc, resultspath, folder)

    doc.build()


def all_results(
        resultspath, 
        filename,
        title = 'Liver-mediated DDI study',
        subtitle = 'All results',
        subject = 'Internal report',
    ):

    print('Creating report..')

    # Cover and title pages
    doc = miblab.Report(
        resultspath,
        filename,
        title = title,
        subtitle = subtitle,
        subject = subject,
    )

    # Two-scan results
    folder = 'twoscan'
    doc.chapter('Two-scan results')
    section_summary(doc, resultspath, folder)
    section_biomarkers(doc, resultspath, folder)
    section_reference(doc, resultspath, folder)
    section_case_notes(doc, resultspath, folder)

    # One-scan results
    folder = 'onescan'
    doc.chapter('One-scan results')
    section_summary(doc, resultspath, folder)
    section_biomarkers(doc, resultspath, folder)
    section_case_notes(doc, resultspath, folder)

    # Secondary results
    doc.chapter('Secondary results')
    section_diurnal(doc, resultspath)
    section_acqtime(doc, resultspath)

    doc.build()


def section_diurnal(doc: miblab.Report, results):

    doc.section('Diurnal variation')

    fig = os.path.join(results, 'twoscan', 'Figures', '_diurnal_function.png')
    caption = (
        "Intra-day changes in hepatocellular uptake (k_he, top row) "
        "and biliary excretion (k_bh, bottom row) of gadoxetate at "
        "for the control (left column) and treatment (right "
        "column). Full lines connect values taken in the same subject "
        "at the same day." 
    )
    doc.figure(fig, caption=caption)



def section_acqtime(doc: miblab.Report, results):

    fig = os.path.join(results, 'onescan_vart', 'Figures', '_effect_plot.png')
    if not os.path.exists(fig):
        return

    doc.section('Acquisition time', clearpage=True)
   
    caption = (
        "The effect of shortening the acquisition time on measured "
        "effect sizes. The figure shows the reference result with 2 "
        "scans, and then repeated results with shorter acquisition "
        "times ranging from 5min to 40min." 
    )
    doc.figure(fig, caption=caption)
            


def section_summary(doc: miblab.Report, results, folder):

    doc.section('Data summary')

    fig = os.path.join(results, folder, 'Figures', '_effect_plot.png')
    caption = (
        "Effect size (%) on hepatocellular uptake (k_he, left) and "
        "biliary excretion (k_bh, right) of gadoxetate. The boxplot "
        "shows median, interquartile range and 95 percent range. The "
        "line plots show individual values for hepatocellular uptake "
        "(k_he, middle) and biliary excretion (k_bh, right) of "
        "gadoxetate of the control (left of plot) and treatment "
        "(right of plot). Grey lines are healthy controls with "
        "rifampicin injection."
    )
    doc.figure(fig, width='7in', caption=caption)

    table = os.path.join(results, folder, 'Analysis', 'k_descriptive_stats.csv')
    caption = (
        "Effect size and absolute values of hepatocellular uptake "
        "(k_he) and biliary excretion (k_bh) of gadoxetate"
    )
    doc.table(table, caption=caption)   
            


def section_reference(doc: miblab.Report, results, folder):

    fig = os.path.join(results, folder, 'Figures', '_compare_to_ref.png')
    if not os.path.exists(fig):
        return

    doc.section('Comparison to reference values', clearpage=True)
    
    caption = (
        "Comparison to reference values in healthy volunteers "
        "treated with the drug and under control conditions."
    )
    doc.figure(fig, caption=caption)

    doc.subsection('Control', clearpage=True)

    table = os.path.join(results, folder, 'Tables', 'reference_liver_control.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "liver biomarkers between this study and healthy "
        "reference data, under control conditions."
    )
    doc.table(table, caption=caption)

    table = os.path.join(results, folder, 'Tables', 'reference_aorta_control.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "aorta biomarkers between this study and healthy "
        "reference data, under control conditions."
    )
    doc.table(table, caption=caption)

    doc.subsection('Treatment', clearpage=True)
        
    table = os.path.join(results, folder, 'Tables', 'reference_liver_drug.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "liver biomarkers between this study and healthy "
        "reference data, under treatment conditions."
    )
    doc.table(table, caption=caption)

    table = os.path.join(results, folder, 'Tables', 'reference_aorta_drug.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "aorta biomarkers between this study and healthy "
        "reference data, under treatment conditions."
    )
    doc.table(table, caption=caption)



def section_biomarkers(doc: miblab.Report, results, folder):

    doc.section('Liver biomarkers', clearpage=True)

    table = os.path.join(results, folder, 'Tables', 'liver_ttest.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "liver biomarkers between control and treatment. The "
        "results are ranked by their p-value, with most significant "
        "differences at the top of the list."
    )
    doc.table(table, caption=caption)

    table = os.path.join(results, folder, 'Tables', 'liver_pairwise.csv')
    caption = (
        "Mean values along with their 95 percent confidence "
        "intervals for all liver biomarkers of the control and "
        "treatment visit. The last column shows the relative change "
        "at the treatment visit. The results are ranked by their "
        "p-value, with most significant differences at the top of the "
        "list."
    )
    doc.table(table, caption=caption)

    doc.section('Systemic biomarkers', clearpage=True)

    table = os.path.join(results, folder, 'Tables', 'aorta_ttest.csv')
    caption = (
        "Results of a pairwise comparison testing for differences in "
        "systemic biomarkers between control and treatment visit. The "
        "results are ranked by their p-value, with most significant "
        "differences at the top of the list."
    )
    doc.table(table, caption=caption)

    table = os.path.join(results, folder, 'Tables', 'aorta_pairwise.csv')
    caption = (
        "Mean values along with their 95 percent confidence "
        "intervals for all systemic biomarkers at the control and "
        "and treatment visit. The last column shows the relative "
        "change at the treatment visit. The results are ranked by "
        "their p-value, with most significant differences at the top "
        "of the list."
    )
    doc.table(table, caption=caption)



def section_case_notes(doc: miblab.Report, results, folder):

    doc.section('Case notes', clearpage=True)

    # Get data
    file = os.path.join(results, folder, 'Analysis', 'parameters_rep.csv')
    subjects = pd.read_csv(file).subject.unique()

    for i, subject in enumerate(subjects):
        if i>0:
            doc.clearpage()
        subj = str(subject).zfill(3)

        doc.subsection('Subject ' + subj)

        # Images
        fig = os.path.join(results, folder, 'Plots', subj +'_control.png')
        caption = (
            "Signal-time curves for subject "+subj+" at the "
            "control visit."
        )
        doc.figure(fig, width='4.5in', caption=caption)

        fig = os.path.join(results, folder,  'Plots', subj +'_drug.png')
        caption = (
            "Signal-time curves for subject "+subj+" at "
            "the treatment visit."
        )
        if os.path.exists(fig):
            doc.figure(fig, width='4.5in', caption=caption)

        # Tables
        doc.clearpage()

        table = os.path.join(results, folder, 'Tables', subj + '_liver.csv')
        caption = "Values for liver of subject "+subj
        doc.table(table, caption=caption)

        table = os.path.join(results, folder, 'Tables', subj + '_aorta.csv')
        caption = "Values for aorta of subject "+subj
        doc.table(table, caption=caption)