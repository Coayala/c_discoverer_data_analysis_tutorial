Tutorial for the analysis of LC-MS/MS feature abundance and annotation
data
================
Christian Ayala

# Introduction

This repository intendts to serve as a guide on how to analyze Feature
abundance and annotation data generated with the **Compound Discoverer**
software. The idea is to expand and customize the analysis and figures
that can be generated directly with Compound Discoverer to allow for
more detailed or in-depth analysis of the data. All the code mentioned
in this tutorial can be found in the file `c_disc_analysis.Rmd` of this
repository.

# 1. Exporting data from Compound Discoverer

In the **Compounds** tab of Compound Discoverer results visualization
select the following *Field* to export:

-   Annot. $\Delta$Mass \[ppm\]
-   Area
-   Calc. MW
-   Formula
-   Gap Status
-   Gap Fill Status
-   MS2
-   Name
-   RT \[min\]

Then, right-click on the headers and select `Export -> As Excel`

<p align="center">
<img src="pics/c_discoverer_fields.png" width="350" title="Selecting Fields in Compound Discoverer">
<img src="pics/c_discoverer_export_excel.png" width="350" title="Export as an Excel File">
</p>

# 2. Analyzing data with RStudio
