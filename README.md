
<!-- README.md is generated from README.Rmd. Please edit that file -->

# FunBootBand

<!-- badges: start -->
<!-- badges: end -->

Compute statistical bands from curve data using a functional approach
and bootstrapping.

At the core, this is an implementation of the method developed by
Sutherland et al. (1988) and Olshen et al. (1989), described in detail
in Lenhoff et al. (1999). The method was originally written as a MATLAB
program by Doris Oriwol and later translated into R and extended with an
approach to handle hierarchical data (see also
<https://github.com/koda86/floa>).

<!-- Bugfix bei der Berechnung -->

More details can be found in this publication (and the vignette):

Koska, D., oriwol, D., & Maiwald, C. (2023). Comparison of statistical
models for characterizing continuous differences between two
biomechanical measurement systems. Journal of Biomechanics, J. Biomech.
149, <https://doi.org/10.1016/j.jbiomech.2023.111506>.

## Citation

You can quote the package …

## Issues

Please report software bugs or other problems by searching existing
issues or creating a new issue here.

## Contributing

If you find issues or bugs feel free to send me an E-Mail
(<daniel.koska@hsw.tu-chemnitz.de>). If you want to fix them yourself,
please do, and submit a pull request so it can be reviewed and merged.

## Licence

## Installation

You can install the development version of FunBootBand from
[GitHub](https://github.com/koda86/FunBootBand) with:

``` r
# install.packages("devtools")
devtools::install_github("koda86/FunBootBand")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(FunBootBand)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub and CRAN.

### References

-   Lenhoff, M.W., Santner, T.J., Otis, J.C., Peterson, M.G., Williams,
    B.J., Backus, S.I., 1999. Bootstrap prediction and confidence bands:
    a superior statistical method for analysis of gait data. Gait
    Posture 9 (1), 10–17. <http://dx.doi.org/10.1016/s0966->
    6362(98)00043-5.

-   Olshen, R.A., Biden, E.N., Wyatt, M.P., Sutherland, D.H., 1989. Gait
    analysis and the bootstrap. Ann. Statist. 17 (4),
    <http://dx.doi.org/10.1214/aos/1176347372>.

-   Sutherland, D., Olshen, R., Biden, E., Wyatt, M., 1988. Development
    of Mature Walking. Mac Keith Press.
