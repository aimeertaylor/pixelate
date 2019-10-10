# pixelate

An R package to pixelate spatial predictions according to the uncertainty that surrounds them.

<!--- e.g. ![ ](link-to-image) --->
<!--- I have not figured out how to export a single code chunck figure from the vignette
and load here since best practice dictates that vignettes should be built on installation;
see https://github.com/r-lib/devtools/issues/584. 
Perhaps we could add some static examples based not on pf. That way we can avoid the need to 
sync and thus avoid any conflict between images here in the manuscript / vugnette? 
---> 

## Prerequisites

The package **pixelate** is an R package. 
It was developed in R version 3.6.1 (2019-07-05) using RStudio. 
To download and install R please go to [cran.r-project.org](https://cran.r-project.org).
To download and install RStudio please go to [rstudio.com](https://rstudio.com/).

## Installation

A development version of **pixelate** is available on Github. 
It can be installed in R using `install_github` from the **devtools** package.
At the time of writing (9th Oct 2019) 
an in-development version of devtools (version 2.2.1.9000)
was needed to build **pixelate**'s vignette (central to understading **pixelate**)
upon installation. Please follow the code below to ensure the vignette builds. Also, 
please ensure `dependencies = TRUE`. If you encounter problems with dependencies, 
see section entitled **Dependency issues** below. 

_**Please be aware**_, installation with `build_vignettes = TRUE` takes several minutes 
(approximately five on a Macbook Pro) 
because the vignette includes plots that are slow to generate. Apologies for the 
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
please install **ggplot2** if required and attach both packages as follows. 

```r
# Install ggplot2 from CRAN if required
if (!require("ggplot2")) install.packages("ggplot2") 

# Attach packages pixelate and ggplot2 
library(pixelate) 
library(ggplot2)
```

Thereafter, please see the documentation of `pixelate` (accessed by `?pixelate` of `help(pixelate)`) 
and read the **pixelate** vignette for quick and detailed examples of **pixelate**'s usage. 
The vignette can be accessed as follows. 

<!--- Avoid examples here s.t. user follows the vignette --->

```r
# Load vignette
vignette("pixelate")
```

## Dependency issues

Originally, we encountered problems with dependencies required for a
single function, `ggplot::coord_sf`, that features in the vignette of **pixelate**. 
Setting `dependencies = TRUE` in `devtools::install_github` appears to solve these problems. 
If problems persist, however, please consider the comments and code below. 

The function `ggplot::coord_sf`, which features in the vignette of **pixelate**
requires the **sf** package that depends on the **units** package. (As an aside, 
please be aware that `ggplot::coord_sf` is not essential: it enables plots to
be enhanced via the addition of shape file details, e.g. country borders). 
The package **sf** is suggested by **ggplot2** and **units** sometimes fails to install. 
If this happens, first try to install **units** as `type = binary` if required. 
If that fails, see instructions for other solutions online 
(e.g. https://community.rstudio.com/t/trouble-installing-packages-in-the-quickstart/23800). 
After installing **units**, install **sf** and **ggplot2** if required. 
If all of the above fail, you can still use `pixelate` and plot its output. 
You just wont be able to add country borders etc. to plots.

```{r setup, eval = FALSE}
# Install packages from CRAN if required (this code is not evaluated)
if (!requireNamespace("units", quietly = TRUE)) install.packages("units", type = 'binary')
if (!requireNamespace("sf", quietly = TRUE)) install.packages("sf")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2", dependencies = TRUE)
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

<!--- Please make sure to update tests as appropriate. --->

## License
[MIT](https://choosealicense.com/licenses/mit/)



<!--- ## Acknowledgements 
Acknowledge everyone who helps test code (e.g. PM)
