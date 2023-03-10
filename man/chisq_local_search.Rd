% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/chisq_tabu_bin_bisect.R
\name{chisq_local_search}
\alias{chisq_local_search}
\title{Finding the Value of Chi Square Statistics for Dynamic Significance Level}
\usage{
chisq_local_search(level_of_significance, dfreedom, tol,
max_iteration, li, ui, is_echo)
}
\arguments{
\item{level_of_significance}{a value between 0 and 1}

\item{dfreedom}{the degree of freedom for given system}

\item{tol}{tolerance quantity, usually 0.001}

\item{max_iteration}{the maximum iteration for the searching}

\item{li}{the lower bound of the search space}

\item{ui}{the upper bound of the search space}

\item{is_echo}{a boolean paramter, by default, is_echo is FALSE}
}
\value{
Return either numeric value or iteration table
}
\description{
\code{chisq_local_search} is used to find the critical value of
chi square statistics at level of significance between 0 and 1
}
\details{
\code{chisq_local_search} function is used to find the
critical value of the chi square statistics for given level of
significance between 0 and 1 with degree of freedom of the system.
This critical value is used to construct the 95\% confidence level
for a certain parameter, and this critical value is used to conduct hypothesis testing as wells.

The value of \code{dfreedom} is positive integers. For example,
in multiple linear regression, the degree of freedom for given
model is n-k-1, \code{tol} is the tolerance quantity, usually, it's
0.001. The tolerance quantity is used to return the critical value
of statistics chi square correct to ith-decimal place.
\code{max_iteration} is the desired number iteration as
termination criteria. \code{li} and \code{ui} are, respectively,
lower bound and the upper bound of the search space, and the
minimum value for \code{li} is 0. \code{is_echo} is a boolean variable, if \code{is_echo} is FALSE, then \code{chisq_local_search} function
return a numeric value only, and if \code{is_echo} is TRUE,
then \code{chisq_local_search} return iteration table.
}
\note{
In order to check either the critical value obtained from \code{chisq_local_search} is, indeed,
the desired critical value, \code{pchisq} function can be used.
}
\examples{
# Given the critical value is 0.025 and degree od freedom is 29,
# then the critical value of the chi square statistics is computed as follows.
chisq_local_search(
level_of_significance = 0.025, dfreedom = 29, tol = 0.001,
max_iteration = 50
) #[1] 45.72261
pchisq(45.72261, df = 29) #[1] 0.9750018
1 - pchisq(45.72261, df = 29) #[1] 0.02499816 ~= 0.025
}
\references{
Howard, J.P. (2017). Computational Methods for Numerical Analysis with R (1st ed.).
Chapman and Hall/CRC. https://doi.org/10.1201/9781315120195

Laguna M. (2018) Tabu Search. In: Mart?? R., Pardalos P., Resende M. (eds) Handbook
of Heuristics. Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_24

Ronald E. Walpole, Raymond H. Myers, Sharon L. Myers, and Keying Ye. Probability
& Statistics for Engineers & Scientists, Nineth Global edition (9th glob. ed.).
Pearson.

Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2009.
Introduction to Algorithms, Third Edition (3rd. ed.). The MIT Press.
}
\seealso{
If one need the critical value of f for dynamic level of
significance, then one can use \code{f_local_search} function.
Similarly, for z and t statistics.
}
