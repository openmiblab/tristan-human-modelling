![example-result](_static/Result.png)

---

# Measuring drug-mediated inhibition of liver transporters

[![Code License: Apache 2.0](https://img.shields.io/badge/License-Apache%202.0-blue.svg?style=flat-square&logo=apache&color=blue)](https://www.apache.org/licenses/LICENSE-2.0) [![Data License: CC BY 4.0](https://img.shields.io/badge/Data%20License-CC%20BY%204.0-lightgrey.svg)](https://creativecommons.org/licenses/by/4.0/) [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15609209.svg)](https://doi.org/10.5281/zenodo.15609209) [![Input Data](https://img.shields.io/badge/input%20data-Zenodo-3776AB?logo=databricks&logoColor=white)](https://doi.org/10.5281/zenodo.15610261) [![Output Data](https://img.shields.io/badge/output%20data-Zenodo-FF8C00?logo=databricks&logoColor=white)](https://doi.org/10.5281/zenodo.15610350)

## 📚 Context 

The liver is responsible for filtering waste products from the blood, and evacuating these by excreting them to bile. Drugs commonly affect these procesesses, potentially leading to toxic side effects.

If a drug inhibits excretion into the bile, then harmful material can get stuck in liver cells and cause damage to the liver. This is commonly referred to as drug-induced liver injury (DILI). If a drug on the other hand blocks the uptake of these materials into the liver, they will circulate in the blood for too long and may cause harm elsewehere.

When dealing with novel drugs, either in basic research or drug development, it is often unknown to what extent a drug affects the uptake or excretion into the liver. This can expose trial participants to significant risk, or even patients in clinical practice if the risk is not identified during development.

In order to mitigate these risks, the [TRISTAN project](https://www.imi-tristan.eu/) developed an MRI-based method to measure the effect of drugs in liver transporters directly. Proof of concept was provided in preclinical studies and pilot studies in humans, showing that the method is able to detect inhibition of uptake and excretion caused by drugs. 

The pipeline in this repository was used to generate the results in humans. It can be used to reproduce them independently, but also to analyse new data acquired in the same way. 

## 🛠️ What does it do? 

The **inputs** to the pipeline are signal-time curves in regions-of-interest in the liver and aorta, using a dynamic gadoxetate-enhanced MRI acquisition, of any number of human subjects. The data must be in .dmr format as in [this repository](https://doi.org/10.5281/zenodo.15610261) with data acquired during the TRISTAN project:

The **outputs** are measurements of gadoxetate uptake rates into hepatocytes, and excretion rates into bile, as well as secondary results such as extracellular volumes or systemic parameters like cardiac output. Key numerical results from the TRISTAN studies are saved in [this database](https://doi.org/10.5281/zenodo.15610350).

Apart from the main biomarkers, the pipeline creates plots, performs statistical analysis and summarises all results in pdf reports. 

## 	📄 Code structure

The **src** folder contains all the source code, with the top level entry scripts, which call on functions in the subfolder **methods**.

The **build** folder contains the output produced by the scripts in **src**. It can be deleted and will be recreated by the scripts. The output of the script **tristan.py** is included in this repository as an example.

## 💻 Usage

The pipeline can be run after installing the requirements:

```console
pip install -r requirements.txt
```

The main scripts in the **src** and all be run independently. They all reproduce results that were generated during the TRISTAN project, except for the script **analyze_newdrug.py** which is a template for future application of the pipeline to newly discovered drugds. The script **analyze_rifampicin.py** only generates the primary results, and the notebook with the same name **analyze_rifampicin.ipynb** is a narrative step-by-step guide to the calculation. 

In order to reproduce existing results, delete the **build** folder and run the script **tristan_all.py**. The entire calculation may take several hours on a laptop computer. The pipeline will download its input data from [a public archive](https://zenodo.org/records/15301607). The local copy of the data is deleted after completion of the analysis. 

In order to analyse a new drug with the same method, perform the following steps:

1. Acquire data as laid out in the [experimental protocol](https://archive.ismrm.org/2024/4015.html) and derive signal-time curves (miblab pipeline for this part coming soon).

2. Save your data in [dmr format](https://openmiblab.github.io/pydmr/format.html) in the same way as the datasets on [the public archive](https://zenodo.org/records/15301607).

3. Then run the script **analyze_newdrug.py** making sure to replace the placeholder values by those that describe your data. When the computation finishes, the results will be added to the build folder.

## ❤️ Citation 

Thazin Min, Marta Tibiletti, Paul Hockings, Aleksandra Galetin, Ebony Gunwhy, Gerry Kenna, Nicola Melillo, Geoff JM Parker, Gunnar Schuetz, Daniel Scotcher, John Waterton, Ian Rowe, and Steven Sourbron. Measurement of liver function with dynamic gadoxetate-enhanced MRI: a validation study in healthy volunteers. [Proc Intl Soc Mag Reson Med, Singapore 2024, 4015](https://archive.ismrm.org/2024/4015.html).

## 💰 Funding 

The work was performed as part of the [TRISTAN project](https://www.imi-tristan.eu/) on imaging biomarkers for drug toxicity. The project was EU-funded through the [Innovative Health Initiative](https://www.ihi.europa.eu/).

[![TRISTAN](_static/tristan-logo.jpg)](https://www.imi-tristan.eu/)

## 👥 Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->
<table>
  <tbody>
    <tr>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/EbonyGunwhy"><img src="https://avatars.githubusercontent.com/u/100941618?v=4" width="100px;" alt="Ebony Gunwhy"/><br /><sub><b>Ebony Gunwhy</b></sub></a><br /></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/EShalom"><img src="https://avatars.githubusercontent.com/u/79933818?v=4" width="100px;" alt="Eve Shalom"/><br /><sub><b>Eve Shalom</b></sub></a><br /></td>
      <td align="center" valign="top" width="14.28%"><a href="https://github.com/plaresmedima"><img src="https://avatars.githubusercontent.com/u/6051075?v=4" width="100px;" alt="Steven Sourbron"/><br /><sub><b>Steven Sourbron</b></sub></a><br /></td>
    </tr>
  </tbody>
</table>

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->

