# pixelate

An R package to pixelate spatial predictions according to the uncertainty that surrounds them.

<!--- e.g. ![ ](link-to-image) --->
<!--- I have not figured out how to export a single code chunck figure from the vignette
and load here since best practice dictates that vignettes should be built on installation;
see https://github.com/r-lib/devtools/issues/584. 
Let's add some static examples based not on pf. That way we can avoid the need to 
sync and thus avoid any conflict between images here in the manuscript / vugnette? 
---> 

## Prerequisites

The package **pixelate** is an R package. It was developed in R version 3.6.1 (2019-07-05) using RStudio. 
To download and install R please go to [cran.r-project.org](https://cran.r-project.org).
To download and install RStudio please go to [rstudio.com](https://rstudio.com/). Please ensure you have the latest version of R or at least version 3.5.0. 

## Installation

A development version of **pixelate** is available on Github. 
It can be installed in R using `install_github` from the **devtools** package.
At the time of writing (9th Oct 2019) 
an in-development version of devtools (version 2.2.1.9000)
was needed to build **pixelate**'s vignette upon installation. 
To ensure **pixelate** installs and the vignette builds, 
please follow the code below and accept any suggested package updates. 

_**Please be aware**_, the vignette takes several minutes to build 
(up to five on an Early 2015 MacVook Pro)
because it includes plots that are slow to generate. Apologies for the 
wait and thank you for your patience. 

<!--- (We chose `build_vignettes = TRUE` over removing `inst/doc` from .gitignore following 
Hadley Wickham's advice; see https://github.com/r-lib/devtools/issues/584). --->

```r
# Step 1) install devtools as required: 

if (!require("devtools")) { # If devtools is not intalled
  
  # Install stable version from CRAN:  
  install.packages("devtools") 
  
  # Extract and compare the stable version: 
  vdetools = as.character(packageVersion("devtools"))
  vcompare = compareVersion(vdetools, '2.2.1.9000')
  
  if (vcompare < 0) {  # If the stable version is < ‘2.2.1.9000’
    
    # Install in-development version from Github
    devtools::install_github("r-lib/devtools") 
  }

} else { # If devtools is already intalled  
  
  # Extract and compare the installed version: 
  vdetools = packageVersion("devtools")
  vcompare = compareVersion(as.character(vdetools), '2.2.1.9000')
  
  if (vcompare < 0) { # If the installed version is < ‘2.2.1.9000’
    
    # Install in-development version from Github
    devtools::install_github("r-lib/devtools") 
  }
}


# Step 2) install pixelate from GitHub 
devtools::install_github("artaylor85/pixelate", build_vignettes = TRUE, dependencies = TRUE)
```

## Usage

The **pixelate** package centres around a single function `pixelate`.
To use `pixelate` and visualise its output following our examples, 
simply load and attach **pixelate** then read the **pixelate** vignette for quick and detailed examples. 
In addition (and if you did not build the vignette upon installation), please see the documentation of `pixelate` 
(accessed by `?pixelate` or `help(pixelate)`).  

```r
# load and attach
library(pixelate) 

# load the vignette for quick and detailed examples
vignette("pixelate") 

# Access documentation for pixelate()
?pixelate 
```

<!--- Avoid examples here s.t. user follows the vignette --->

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

<!--- Please make sure to update tests as appropriate. --->

## License
[MIT](https://choosealicense.com/licenses/mit/)


## Acknowledgements 
Thank you to Pierre Jacob, Pamela Martinez Vargas, Rene Niehu and Pablo Martínez de Salazar for help testing package installation. 
