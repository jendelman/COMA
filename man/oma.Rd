% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/oma.R
\name{oma}
\alias{oma}
\title{Optimal Mate Allocation}
\usage{
oma(
  dF,
  parents,
  matings,
  ploidy,
  K,
  tol = 1e-06,
  dF.adapt = NULL,
  solver = "ECOS"
)
}
\arguments{
\item{dF}{inbreeding rate}

\item{parents}{parents data frame (see Details)}

\item{matings}{matings data frame (see Details)}

\item{ploidy}{ploidy}

\item{K}{kinship matrix}

\item{tol}{tolerance, values below this set to 0}

\item{dF.adapt}{see Details}

\item{solver}{solver for CVXR (default is "ECOS")}
}
\value{
list containing
\describe{
\item{response}{data.frame with realized dF, merit, n.parent, n.mate}
\item{oc}{data frame of optimal contributions for each individual}
\item{om}{data frame of optimal allocations for each mating}
}
}
\description{
Optimize the allocation for each mating
}
\details{
The first three columns of \code{parents} should be named "id", "min", "max", with an optional fourth column "female" to indicate sex in dioecious species. Additional columns can be used to specify constraints on linear combinations of the contributions. The values are the coefficients of the contribution variables, and the name of each column specifies the right-hand side of the constraint. Each name must begin with "lt","gt", or "eq", followed by a non-negative numeric value. For example, "lt0.5" means less than or equal to 0.5.

The data.frame \code{matings} has five columns: "female, male, merit, min, max" for dioecious species, or else "parent1, parent2, merit, min, max".

The average inbreeding coefficient of the current generation is based on all individuals in \code{K}, which may exceed the list of individuals in \code{parents}.

It is possible that no feasible solution exists for the specified \code{dF}. Argument \code{dF.adapt} can be used to automatically progressively higher values.The software increases the dF limit by dF.adapt$step up to the smaller of dF.adapt$max or the realized value under the original dF, in an attempt to find a solution with less inbreeding.
}
