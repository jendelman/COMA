R/COMA
================
Jeffrey Endelman

Breeders have long appreciated the need to balance selection for
short-term genetic gain with maintaining genetic variance for long-term
gain. The COMA package implements selection strategies known as Optimum
Contribution Selection (OCS) and Optimum Mate Allocation (OMA). OCS
maximizes the average genomic-estimated breeding value (GEBV) of the
parents, weighted by their contribution to the next generation. OMA
maximizes the average genomic prediction of mate performance (GPMP),
weighted by the contribution of each mating, which is called mate
allocation. Constraints on inbreeding rate are used to ensure genetic
variance is not depleted too quickly. The “C” in COMA stands for Convex
because the software exploits the convex nature of the problem to
efficiently find the global optimum.

Financial support for developing COMA has come from the USDA National
Institute of Food and Agriculture (NIFA) Award 2020-51181-32156. Please
cite the [manuscript](https://doi.org/10.1093/genetics/iyae193) if you
use it. [Vignette 1](https://jendelman.github.io/COMA/Vignette1.html)
provides examples of using the software with data from the University of
Wisconsin-Madison potato breeding program.

To install and load the package:

``` r
install.packages("devtools")
devtools::install_github("jendelman/COMA", build_vignettes=FALSE)
library(COMA)
```

For a complete specification of package functions, consult the
[reference manual.](https://jendelman.github.io/COMA/manual.pdf)
