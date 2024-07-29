# LocationSmoothedBART

## Overview

LocationSmoothedBART (lsBART) is an innovative R package that extends the Bayesian Additive Regression Trees (BART) model. It is specifically designed for scalar on function regression and is particularly suited for medical image analysis in radiation treatment planning. This package emphasizes identifying significant regions within predictor functions, ensuring accurate and reliable contour delineation.

## Features

- Extends the traditional BART model for more precise medical image analysis.
- Utilizes shape features for contour quality assurance, independent of the imaging platform.
- Incorporates Shapley values for detailed error analysis in contour delineation.

## Installation

``` r
# To install the latest version from GitHub:
devtools::install_github("wootz101/LocationSmoothedBART")
```

## Documentation

Vignette for Continuous and Binary response variable is available at https://wootz101.github.io/LocationSmoothedBART/

For details on functions, please call the help function with package name.

``` r
help("LocationSmoothedBART")
```

## Contributing

Contributions to the LocationSmoothedBART package are welcome. Please refer to our contributing guidelines for more information.

## Citation

If you use LocationSmoothedBART in your research, please cite:

Wooten, Z. et al. (2024). "Location smoothed BART: A method for interpretable and robust quality assurance of organ contours in radiotherapy treatment planning" (Under Review)

This work is built on the BART R package by Sparapani, Spanbauer, and McCulloch. For additional background information and related methodologies, refer to the BART R package documentation and associated research papers.

  Sparapani R, Spanbauer C, McCulloch R (2021). “Nonparametric Machine Learning and Efficient Computation with Bayesian
  Additive Regression Trees: The BART R Package.” _Journal of Statistical Software_, *97*(1), 1-66.
  doi:10.18637/jss.v097.i01 <https://doi.org/10.18637/jss.v097.i01>.

For the implementing the TreeShap algorithm within the BART framework we built on the TreeShap package by Komisarczyk, Kozminski, Maksymiuk, and Biecek. For additional background information and related methodologies, refer to the TreeShap package documentation and associated research papers.

  Komisarczyk K, Kozminski P, Maksymiuk S, Biecek P (2023). _treeshap: Fast SHAP values computation for tree ensemble
  models_. R package version 0.1.1.9001, <https://github.com/ModelOriented/treeshap>.


## Contact

Zachary Wooten
ztw5@rice.edu
(https://wootz101.github.io/Website2/)
