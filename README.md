# QMC with one categorical variable

This repository contains the R code used to produce the results for the two examples presented in our paper ["QMC with one categorical variable"](https://www.arxiv.org/abs/2506.16582).

## Repository Structure

- `Toy_Example/`  
  Contains files for the toy example:
  - `toy_driver.R` – Main script to run the toy example.  
  - `toy_functions.R` – Functions used to compute the proposed estimator using MC and RQMC points.  
  - `plot_results_toy.R` – Functions for plotting the variance curves shown in the paper.  

- `Flood_example/`  
  Contains files for the Saint-Venant flood example:
  - `flood_driver.R` – Main script to run the example.  
  - `flood_model.R` – Functions used to compute the proposed estimator using MC and RQMC points.  
  - `plot_results_flood.R` – Functions for plotting the variance curves shown in the paper.  

- `rsobol.R`, `sobol_Cs.col`, `fiftysobol.col`  
 Shared functions used by both examples to generate (scrambled) Sobol point sets, from Art Owen's original implementation (https://artowen.su.domains/code/).

## Getting Started

1. Clone the repository:

```bash
git clone https://github.com/yourusername/qmc4categorical.git
```

2. Open R or RStudio and set the working directory to the repository root:

```r
setwd("path/to/qmc4categorical")
```

3. Install ggplot2 and VGAM (if not already installed)

```r
install.packages("ggplot2")
install.packages("VGAM")
```

## Running the examples

- To run the toy example, run in R or RStudio:
```r
source("Toy_Example/toy_driver.R")
```

- To run the flood model example, run in R or RStudio:
```r
source("Flood_Example/flood_driver.R")
```

The plots will be saved in the corresponding example folders (as `toy_results` and `flood_results` folders).

