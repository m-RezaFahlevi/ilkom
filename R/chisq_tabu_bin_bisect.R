#' Finding the Value of Chi Square Statistics for Dynamic Significance Level
#'
#' \code{chisq_local_search} is used to find the critical value of
#' chi square statistics at level of significance between 0 and 1
#'
#' @usage chisq_local_search(level_of_significance, dfreedom, tol,
#' max_iteration, li, ui, is_echo)
#'
#' @param level_of_significance a value between 0 and 1
#' @param dfreedom the degree of freedom for given system
#' @param tol tolerance quantity, usually 0.001
#' @param max_iteration the maximum iteration for the searching
#' @param li the lower bound of the search space
#' @param ui the upper bound of the search space
#' @param is_echo a boolean paramter, by default, is_echo is FALSE
#'
#' @details \code{chisq_local_search} function is used to find the
#' critical value of the chi square statistics for given level of
#' significance between 0 and 1 with degree of freedom of the system.
#' This critical value is used to construct the 95\% confidence level
#' for a certain parameter, and this critical value is used to conduct hypothesis testing as wells.
#'
#' The value of \code{dfreedom} is positive integers. For example,
#' in multiple linear regression, the degree of freedom for given
#' model is n-k-1, \code{tol} is the tolerance quantity, usually, it's
#' 0.001. The tolerance quantity is used to return the critical value
#' of statistics chi square correct to ith-decimal place.
#' \code{max_iteration} is the desired number iteration as
#' termination criteria. \code{li} and \code{ui} are, respectively,
#' lower bound and the upper bound of the search space, and the
#' minimum value for \code{li} is 0. \code{is_echo} is a boolean variable, if \code{is_echo} is FALSE, then \code{chisq_local_search} function
#' return a numeric value only, and if \code{is_echo} is TRUE,
#' then \code{chisq_local_search} return iteration table.
#'
#' @return Return either numeric value or iteration table
#'
#' @note In order to check either the critical value obtained from \code{chisq_local_search} is, indeed,
#' the desired critical value, \code{pchisq} function can be used.
#'
#' @references Howard, J.P. (2017). Computational Methods for Numerical Analysis with R (1st ed.).
#' Chapman and Hall/CRC. https://doi.org/10.1201/9781315120195
#'
#' Laguna M. (2018) Tabu Search. In: Martí R., Pardalos P., Resende M. (eds) Handbook
#' of Heuristics. Springer, Cham. https://doi.org/10.1007/978-3-319-07124-4_24
#'
#' Ronald E. Walpole, Raymond H. Myers, Sharon L. Myers, and Keying Ye. Probability
#' & Statistics for Engineers & Scientists, Nineth Global edition (9th glob. ed.).
#' Pearson.
#'
#' Thomas H. Cormen, Charles E. Leiserson, Ronald L. Rivest, and Clifford Stein. 2009.
#' Introduction to Algorithms, Third Edition (3rd. ed.). The MIT Press.
#'
#' @seealso If one need the critical value of f for dynamic level of
#' significance, then one can use \code{f_local_search} function.
#' Similarly, for z and t statistics.
#'
#' @importFrom stats pchisq runif
#'
#' @examples
#' # Given the critical value is 0.025 and degree od freedom is 29,
#' # then the critical value of the chi square statistics is computed as follows.
#' chisq_local_search(
#' level_of_significance = 0.025, dfreedom = 29, tol = 0.001,
#' max_iteration = 50
#' ) #[1] 45.72261
#' pchisq(45.72261, df = 29) #[1] 0.9750018
#' 1 - pchisq(45.72261, df = 29) #[1] 0.02499816 ~= 0.025
#' @export
chisq_local_search <- function(level_of_significance, dfreedom, tol, max_iteration, li = 0, ui = 50, is_echo = FALSE) {
    given_area <- 1 - level_of_significance # a.k.a confidence level

    # data gathering
    li_vect <- c()
    ui_vect <- c()
    chisq_val_vect <- c()
    areavect <- c()
    tol_vect <- c()

    iter <- 0
    current_tol <- abs(li - ui)
    while (current_tol > tol && iter <= max_iteration) {
        tol_vect <- append(tol_vect, current_tol) # gather the tolerance quantity

        iter <- iter + 1
        chisq_val_estimate <- runif(n = 1, min = li, max = ui)
        calculated_area <- pchisq(q = chisq_val_estimate, df = dfreedom)

        # data gathering
        li_vect <- append(li_vect, li)
        ui_vect <- append(ui_vect, ui)
        chisq_val_vect <- append(chisq_val_vect, chisq_val_estimate)
        areavect <- append(areavect, calculated_area)

        # update the tabu by reducing the interval of
        # pseudo random number generator
        ifelse(
            calculated_area > given_area,
            ui <- chisq_val_estimate, li <- chisq_val_estimate
        )
        current_tol <- abs(li - ui)
    }

    chisq_interpolate <- (li + ui) / 2 # interpolate the chisquare value

    # data gathering
    li_vect <- append(li_vect, li)
    ui_vect <- append(ui_vect, ui)
    chisq_val_vect <- append(chisq_val_vect, chisq_interpolate)
    areavect <- append(areavect, pchisq(q = chisq_interpolate, df = dfreedom))
    tol_vect <- append(tol_vect, abs(li - ui))

    # create iteration table as return value if is_echo is TRUE
    iteration_table <- data.frame(
        "iteration" = seq_len(length(chisq_val_vect)),
        "lower_bound" = li_vect,
        "upper_bound" = ui_vect,
        "chisq_estimate" = chisq_val_vect,
        "cumulative_probability" = areavect,
        "absolute_error" = tol_vect
    )

    # return phase. Is is_echo is false, then return only chisq_value
    ifelse(
        is_echo,
        return(iteration_table), return(chisq_val_vect[length(chisq_val_vect)])
    )
}
