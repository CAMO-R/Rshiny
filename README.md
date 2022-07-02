# CAMO
CAMO is an analytical software with R Shiny based graphical user interface (GUI) for evaluating omics congruence of model organisms. It performs threshold-free Bayesian differential analysis and generates quantitative concordance and discordance scores (c-scores and d-scores) both genome-wide and at pathway level for all pair-wise studies. Based on the c-scores/d-scores, CAMO conducts a series of downstream machine learning and bioinformatics analysis with interactive visualization for pathway knowledge retrieval and topological gene module detection.

## How to start the app

#### Requirement
* R >= 4.0.0
* Rcpp >= 1.0.0
* Shiny >= 1.0.0

MacOS users may need to install xcode first. Run the following code in terminal to install:
```
xcode-select --install
```
MacOS users may encounter error ``*** C++11 compiler required; enable C++11 mode in your compiler, or use an earlier version of Armadillo", which can be solved by run following code command in R Console:
```
Sys.setenv("PKG\_CXXFLAGS"="-std=c++11")
```

#### Install the Shiny software
1. Install the CAMO R package and dependency packages following the instruction at [https://github.com/CAMO-R/Rpackage](https://github.com/CAMO-R/Rpackage).
2. Download the CAMO Shiny project at [https://github.com/CAMO-R/Rshiny](https://github.com/CAMO-R/Rshiny) by clicking on "code > Download ZIP" and extract to a local folder renamed as "Rshiny".

#### Start the Shiny software
1. Open "RunShiny.R" file in R console.
2. Set the working directory of R to the directory "path\_to\_Rshiny/" which contains the Shiny project folder "Rshiny" saved above.
3. Click on the "Run App" button.
Note that the installation progress of R packages may take up to a few minutes. Please check the progress in R console and may need to select whether to update all/some/none packages. After all packages has been installed, the CAMO Shiny app will automatically open in your default browser.

## Where to find the full tutorial 
After starting CAMO, users need to select a local directory as the working directory in the Setting page before any analysis. Please refer to the full tutorial for details at 
[https://github.com/CAMO-R/other/blob/main/Rshiny_tutorial/CAMO_RShiny_Tutorial.pdf](https://github.com/CAMO-R/other/blob/main/Rshiny_tutorial/CAMO_RShiny_Tutorial.pdf)


