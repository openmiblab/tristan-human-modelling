import pylatex as pl
from pylatex.utils import NoEscape
from miblab import report

import tristan.report as rep


def build(filename, resultspath):

    print('Creating report..')

    # Cover and title pages
    doc = report.setup(
        resultspath,
        title = 'Two-compound study',
        subtitle = filename + ' results',
        subject = 'D2.10 - Internal report',
    )
    doc.append(pl.NewPage())
    doc.append(NoEscape('\\tableofcontents'))
    doc.append(NoEscape('\\mainmatter'))

    # Two-scan results
    folder = 'twoscan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Two-scan results'))
    rep.section_summary(doc, resultspath, folder)
    rep.section_biomarkers(doc, resultspath, folder)
    rep.section_reference(doc, resultspath, folder)
    rep.section_case_notes(doc, resultspath, folder)

    # One-scan results
    folder = 'onescan'
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'One-scan results'))
    rep.section_summary(doc, resultspath, folder)
    rep.section_biomarkers(doc, resultspath, folder)
    rep.section_case_notes(doc, resultspath, folder)

    # Secondary results
    doc.append(NoEscape('\\clearpage'))
    doc.append(pl.Command('chapter', 'Secondary results'))
    rep.section_diurnal(doc, resultspath)
    #section_acqtime(doc, results)

    report.build(doc, 'report_' + filename, resultspath)

