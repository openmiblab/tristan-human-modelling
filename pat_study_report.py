import pylatex as pl
from pylatex.utils import NoEscape
from miblab import report

import tristan.report as rep


def build(results):

    filename = 'gothenburg'

    print('Creating report..')

    # Cover and title pages
    doc = report.setup(
        results,
        title = 'Gothenburg patient study',
        subtitle = 'Interim analysis',
        subject = 'D2.07 - Internal report',
    )
    doc.append(pl.NewPage())
    doc.append(NoEscape('\\tableofcontents'))
    doc.append(NoEscape('\\mainmatter'))

    # Two-scan results
    folder = 'twoscan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Two-scan results'))
    rep.section_summary(doc, results, folder)
    rep.section_biomarkers(doc, results, folder)
    rep.section_reference(doc, results, folder)
    rep.section_case_notes(doc, results, folder)

    # One-scan results
    folder = 'onescan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'One-scan results'))
    rep.section_summary(doc, results, folder)
    rep.section_biomarkers(doc, results, folder)
    rep.section_case_notes(doc, results, folder)

    # Secondary results
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Secondary results'))
    rep.section_diurnal(doc, results)
    rep.section_acqtime(doc, results)

    report.build(doc, 'report_' + filename, results)


