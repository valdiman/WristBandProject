## README file
# Silicone Wristbands as Personal Passive Samplers for Airborne PCBs

The scripts included here are part of the study “Advances in Understanding Silicone Wristbands
as Tools for Estimating Personal Exposure to Airborne PCB Concentrations” and were developed
to organize and visualize data, estimate sampling rates for static, dynamic, and worn silicone
wristbands (WBs), and perform statistical analyses, including t-tests and similarity tests such
as cosine theta, supporting the use of WBs as personal passive samplers for airborne PCBs.


----------------------
General Information
----------------------

Deposit Title: Silicone Wristbands as Personal Passive Samplers for Airborne PCBs

Contributor information:

Andres Martinez, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

This README file was generated on April 28 2025 by Andres Martinez.

This work was supported by the National Institutes of Environmental Health Sciences (NIH P42ES013661)
and the University of Iowa Environmental Health Sciences Research Center (NIH P30ES005605).
The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data;
in creation of the dataset; and/or in the decision to submit this data for publication or deposit it in a repository.

Subjects: polychlorinated biphenyls; airborne; sampling rate; GC-MS/MS; personal passive sampling


--------
PREREQUISITES & DEPENDENCIES
--------

This section of the ReadMe file lists the necessary software required to run codes in "R".

Software:

Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file: https://www.rstudio.com/products/rstudio/download/

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".
"R Studio" is an open source (no product license required) integrated development environment (IDE) for
"R" and completely free to use. To install "R Studio" follow the instructions below:

Visit the following web address: https://www.rstudio.com/products/rstudio/download/
Click the "download" button beneath RStudio Desktop
Click the button beneath "Download RStudio Desktop". This will download the correct installation file based on the operating system detected.
Run the installation file and follow on-screen instructions.


--------
R FILES AND STRUCTURE
--------

It is recommended to create a project in R (e.g., WristbandProject.Rproj).
Download the project file (.Rproj) and the R subfolder where the scripts are located,
and the Subfolders.R file. Run first the Subfolder.R file, which will generate all the
subfolders for this project. The structure of this project includes an R subfolder where
all the R scripts are located, as previously indicated. There is a Data subfolder where
the data are storage ("Data/IRO"), and then an Output subfolder, where the results are located.


--------
Data
--------

The data files need to run the scripts can be found and downloaded at:
Martinez, A., Metwali, N., Arp, L., Moran, M. E., Marek, R. F., Thorne, P. S. (2026). Dataset for: Individual PCB Measurements in Silicone Wristbands
Used to Estimate Personal Exposure to Airborne PCB Concentrations [Dataset]. University of Iowa. https://doi.org/10.25820/data.007569

The files are:

03_SampleWBMassStudy1.csv

05_SamplePUFConcStudy1.csv

07_SampleWBMassStudy2.csv

08_BlankWBMassStudy3_4_5.csv

09_SampleWBMassStudy3_4_5.csv

10_FrederiksenWBMass2022.csv

11_FrederiksenConc2022.csv

12_logKoa.csv





