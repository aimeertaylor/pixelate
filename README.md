# pixelate

An R package to pixelate spatial predictions according to the uncertainty that surrounds them.

<!--- e.g. ![ ](link-to-image) --->
<!--- Thus far, I have not figured out how to export a single code chunck figure into vignette 
<!--- to load here ---> 
<!--- The above figure shows pixelated predictions of *Plasmodium falciparum* incidence in 2017 in sub-Saharan Africa.
The unpixelated predictions and uncertainty measures were downloaded from the Malaria Atlas Project (www.map.ox.ac.uk). 
They feature in the paper @weiss2019. --->

## Prerequisites

The package **pixelate** is an R package. 
It was developed in R version 3.6.1 (2019-07-05) using RStudio. 
To download and install R please go to [cran.r-project.org](https://cran.r-project.org).
To download and install RStudio please go to [rstudio.com](https://rstudio.com/).

## Installation

A development version of **pixelate** is available on Github 
and can be installed in R using `install_github` from the **devtools** package as follows.

```r
# Install devtools from CRAN if required
if (!require("devtools")) install.packages("devtools") 
 
# Install pixelate from GitHub 
devtools::install_github("artaylor85/pixelate")
```

## Usage

The **pixelate** package centres around a single function `pixelate` 
whose output we visualise using functions from the **ggplot2** package. 
To use `pixelate` and visualise its output following our examples, 
please install ggplot2 if required and attach both packages as follows. 

```r
# Install ggplot2 from CRAN if required
if (!require("ggplot2")) install.packages("ggplot2") 

# Attach packages pixelate and ggplot2 
library(pixelate) 
library(ggplot2)
```

Thereafter, please read the **pixelate** vignette for quick and detailed examples of **pixelate**'s usage. 
The vignette can be accessed as follows. 
Please be aware, it takes a few moments to build.  

<!--- To-add: the figures in the manuscript are fully reproducible following the detailed example. --->
<!--- Avoid examples here s.t. users read following the vignette --->

```r
# Load vignette
vignette("pixelate")
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

<!--- Please make sure to update tests as appropriate. --->

## License
[MIT](https://choosealicense.com/licenses/mit/)

<!--- ## Acknowledgements 
Acknowledge everyone who helps test code (e.g. PM)
