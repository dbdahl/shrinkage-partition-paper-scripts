<!--HOW TO COMPLETE THIS FORM:-->

<!--
1. Checkboxes in this document appear as follows: 

- [ ] This is a checkbox 

To check a checkbox, replace [ ] by [x], as follows: 

- [x] This is a checked checkbox 

Note that older versions of RStudio (versions lower than 1.3) may not create a formatted checkbox but will leave the original characters, i.e., literally "[ ]" or "[x]". It's fine to submit a PDF in this form.
 
2. For text answers, simply type the relevant text in the areas indicated. A blank line starts a new paragraph. 
 
3. Comments (like these instructions) provide additional instructions throughout the form. There is no need to remove them; they will not appear in the compiled document. 

4. If you are comfortable with Markdown syntax, you may choose to include any Markdown-compliant formatting in the form. For example, you may wish to include R code chunks and compile this document in R Markdown.
-->

This form documents the artifacts associated with the article (i.e., the data and code supporting the computational findings) and describes how to reproduce the findings.


# Part 1: Data

- [ ] This paper does not involve analysis of external data (i.e., no data are used or the only data are generated by the authors via simulation in their code).

<!--
If box above is checked and if no simulated/synthetic data files are provided by the authors, please skip directly to the Code section. Otherwise, continue.
-->

- [x] I certify that the author(s) of the manuscript have legitimate access to and permission to use the data used in this manuscript.

<!-- If data are simulated using random number generation, please be sure to set the random number seed in the code you provide -->

## Abstract

<!--
Provide a short (< 100 words), high-level description of the data
-->

The data were obtained from [IPUMS CPS](https://cps.ipums.org/cps/). In March
of each year, a cross-section of individuals from the 50 United States and
the District of Columbia (D.C.) are surveyed regarding their employment and
demographic  information.  Data with harmonized variables are available from the
years 1994 to 2020.  Applying sensible filters for our application (e.g., only
including working individuals) yields 139,555 observations.  For simplicity,
sampling weights are ignored for this  demonstration.

## Availability


- [x] Data **are** publicly available.
- [ ] Data **cannot be made** publicly available.

If the data are publicly available, see the *Publicly available data* section. Otherwise, see the *Non-publicly available data* section, below.

### Publicly available data

- [ ] Data are available online at: https://github.com/dbdahl/shrinkage-partition-paper-scripts
      See in particular the script '0010-clean.R' to read the data from the 'data-raw' directory.

- [ ] Data are available as part of the paper’s supplementary material.

- [ ] Data are publicly available by request, following the process described here:

- [ ] Data are or will be made available through some other mechanism, described here:


<!-- If data are available by request to the authors or some other data owner, please make sure to explain the process of requesting access to the data. -->

### Non-publicly available data

<!--
The Journal of the American Statistical Association requires authors to make data accompanying their papers available to the scientific community except in cases where: 1) public sharing of data would be impossible, 2) suitable synthetic data are provided which allow the main analyses to be replicated (recognizing that results may differ from the "real" data analyses), and 3) the scientific value of the results and methods outweigh the lack of reproducibility.

Please discuss the lack of publicly available data. For example:
-	why data sharing is not possible,
-	what synthetic data are provided, and 
-	why the value of the paper's scientific contribution outweighs the lack of reproducibility.
-->

## Description

### File format(s)

<!--
Check all that apply
-->
- [ ] CSV or other plain text.
- [ ] Software-specific binary format (.Rda, Python pickle, etc.): pkcle
- [ ] Standardized binary format (e.g., netCDF, HDF5, etc.): 
- [x] Other (please specify): The data is in IMPUS's format and read into R using their CRAN package "ipumsr" as shown in the '0010-clean.R' script. 

### Data dictionary

<!--
A data dictionary provides information that allows users to understand the meaning, format, and use of the data.
-->

- [ ] Provided by authors in the following file(s):
- [x] Data file(s) is(are) self-describing (e.g., netCDF files)
- [ ] Available at the following URL: 

### Additional Information (optional)

<!-- 
OPTIONAL: Provide any additional details that would be helpful in understanding the data. If relevant, please provide unique identifier/DOI/version information and/or license/terms of use.
-->

The IMPUS CPS license agreement says, "You may publish a subset of the data to meet journal requirements for accessing data related to a particular publication. Contact us for permission for any other redistribution; we will consider requests for free and commercial redistribution."

# Part 2: Code

## Abstract

<!--
Provide a short (< 100 words), high-level description of the code. If necessary, more details can be provided in files that accompany the code. If no code is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

The code consi

## Description

### Code format(s)

<!--
Check all that apply
-->
- [x] Script files
    - [x] R
    - [ ] Python
    - [ ] Matlab
    - [ ] Other: 
- [x] Package
    - [x] R
    - [ ] Python
    - [ ] MATLAB toolbox
    - [ ] Other: 
- [ ] Reproducible report 
    - [ ] R Markdown
    - [ ] Jupyter notebook
    - [ ] Other:
- [x] Shell script
- [ ] Other (please specify): 

### Supporting software requirements

<!--
Please cite all software packages in the References Section in similar fashion to paper citations, citing packages that are foundational to the research outcome (including packages that implement methods to which you compare your methods). You may elect to not cite packages used for supporting purposes. For R packages, note that running `citation('name_of_package')` often shows how the package authors wish to be cited. 
-->

Software implementing our SP distribution is available as an
R package based on Rust (\url{https://github.com/dbdahl/gourd-package}).

#### Version of primary software used

<!--
(e.g., R version 3.6.0)
-->

R version 4.4.1

#### Libraries and dependencies used by the code

<!--
Include version numbers (e.g., version numbers for any R or Python packages used)
-->

The follow R packages are used in the replication scripts:
+ mvtnorm (1.2.5)
+ coda (0.19-4.1)
+ ipumsr (0.8.1)
+ salso (0.3.38)
+ MASS (7.3-61)
+ fields (16.2)
+ txtplot (1.0.4)

### Supporting system/hardware requirements (optional)

<!--
OPTIONAL: System/hardware requirements including operating system with version number, access to cluster, GPUs, etc.
-->

Reproducing all the results requires running many different R scripts with
different input parameters. The jobs to run by R are generated by shell scripts
written in [fish](https://fishshell.com/) version 3.7.1. In all, there are many
thousands of single-thread jobs to run.  Practically speaking, to reproduce
every result in the paper, one or more Linux servers running for several days
would be need.

### Parallelization used

- [x] No parallel code used
- [ ] Multi-core parallelization on a single machine/node
    - Number of cores used: 
- [ ] Multi-machine/multi-node parallelization 
    - Number of nodes and cores used: 

### License

- [x] MIT License (default)
- [ ] BSD 
- [ ] GPL v3.0
- [ ] Creative Commons
- [ ] Other: (please specify)


### Additional information (optional)

<!--
OPTIONAL: By default, submitted code will be published on the JASA GitHub repository (http://github.com/JASA-ACS) as well as in the supplementary material. Authors are encouraged to also make their code available in a public code repository, such as on GitHub, GitLab, or BitBucket. If relevant, please provide unique identifier/DOI/version information (e.g., a Git commit ID, branch, release, or tag). If the code and workflow are provided together, this section may be omitted, with information provided in the "Location" section below.
-->

The code and workflow are provided together.  See the "Location" section below.

# Part 3: Reproducibility workflow

<!--
The materials provided should provide a straightforward way for reviewers and readers to reproduce analyses with as few steps as possible. 
-->

## Scope

The provided workflow reproduces:

- [x] Any numbers provided in text in the paper
- [ ] The computational method(s) presented in the paper (i.e., code is provided that implements the method(s))
- [x] All tables and figures in the paper
- [ ] Selected tables and figures in the paper, as explained and justified below:

## Workflow

### Location

The workflow is available:

<!--
Check all that apply, and in the case of a Git repository include unique identifier, such as specific commit ID, branch, release, or tag.
-->
- [ ] As part of the paper’s supplementary material.
- [x] In this Git repository:  https://github.com/dbdahl/shrinkage-partition-paper-scripts
- [ ] Other (please specify):

<!--
Indicate where the materials (generally including the code, unless in a separate location and indicated in the previous section) are available. We strongly encourage authors to place their materials (but not large datasets) in a Git repository hosted on a site such as GitHub, GitLab, or BitBucket. If the repository is private during the review process, please indicate the location where it will be available publicly upon publication, and also include the materials as a zip file (e.g., obtained directly from the Git hosting site) as supplementary materials.
-->


### Format(s)

<!--
Check all that apply
-->
- [ ] Single master code file 
- [x] Wrapper (shell) script(s)
- [ ] Self-contained R Markdown file, Jupyter notebook, or other literate programming approach
- [ ] Text file (e.g., a readme-style file) that documents workflow
- [ ] Makefile
- [ ] Other (more detail in *Instructions* below)

### Instructions

<!--
Describe how to use the materials provided to reproduce analyses in the manuscript. Additional details can be provided in file(s) accompanying the reproducibility materials. If no workflow is provided, please state this and say why (e.g., if the paper contains no computational work).
-->

Steps:
+ Clone the Git repository at https://github.com/dbdahl/shrinkage-partition-paper-scripts.
+ Sort the scripts whose file names start with 4 digits, e.g. '0010-clean.R', '0014-jobs-ols', ..., '0230-summarize-simulation.R'
+ Run each script, in numerical order.
    + Scripts whose name end in '.R' should be run by R using, e.g., 'R CMD BATCH 0010-clean.R'.
    + Scripts whose name match the glob '*-jobs-*' produce commands, one per line, that need to be executed before running the next script.
        + You may wish to run these concurrently using, e.g. GNU Parallel (https://www.gnu.org/software/parallel/).
+ In all, there are many thousands of single-thread jobs to run.
    + Practically speaking, to reproduce every result in the paper, many Linux servers running for several days would be need.
    + To avoid needing to run a particular script, the reader may refer to the associated '*.Rout' files that are available in the repository.

### Expected run-time

Approximate time needed to reproduce the analyses on a standard desktop machine:

- [ ] < 1 minute
- [ ] 1-10 minutes
- [ ] 10-60 minutes
- [ ] 1-8 hours
- [ ] > 8 hours
- [x] Not feasible to run on a desktop machine, as described here:

Reproducing all the results requires running many different R scripts with different input parameters.
The jobs to run by R are generated by shell scripts written in [fish](https://fishshell.com/) version 3.7.1.
In all, there are many thousands of single-thread jobs to run.  Practically speaking, to reproduce every result in the paper, many Linux servers running for several days would be need.

### Additional information (optional)

<!--
OPTIONAL: Additional documentation provided (e.g., R package vignettes, demos or other examples) that show how to use the provided code/software in other settings.
-->

# Notes (optional)

<!--
OPTIONAL: Any other relevant information not covered on this form. If reproducibility materials are not publicly available at the time of submission, please provide information here on how the reviewers can view the materials.
-->
