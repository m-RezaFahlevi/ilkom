#' Finding the Value of f Statistics for Dynamics Level of Significance
#'
#' \code{f_local_search} is used to find the critical value of the
#' f statistics for given any level of significance between 0 and 1
#'
#' @usage f_local_search(level_of_significance, dfreedom1, dfreedom2,
#' tol, max_iteration, li, ui, is_echo)
#'
#' @param level_of_significance a value between 0 and 1
#' @param dfreedom1 degree of freedom for the first kind
#' @param dfreedom2 degree of freedom for the second kind
#' @param tol tolerance quantity
#' @param max_iteration the desired number of iteration
#' @param li the lower bound of the search space, by default, it's set to 0
#' @param ui the upper bound of the search space, by default, it's set to 50
#' @param is_echo a boolean parameter, by default, it's FALSE
#'
#' @details \code{f_local_search} function is used to find the
#' critical value of the f statistics for given any level of
#' significance between 0 and 1. This critical value is used to
#' construct the 95\% confidence interval for certain parameters,
#' and this critical value is used in hypothesis testing as wells.
#' \code{f_local_searc} function contain 8 parameters with the
#' following details, 3 parameters, which by defaults, already
#' contain a value, and 5 parameters are not.
#'
#' \code{level_of_significance} is a value between 0 and 1, and it's
#' corresponds to confidence value, if confidence limits is known,
#' then the level of significance can be obtained by subtracting
#' confidence limits from 1. \code{dfreedom1} is the degree of
#' freedom of the first kind, similarly for \code{dfreedom2}. For
#' example, in the hypothesis testing concerning either certain
#' parameters in the multiple linear regression model is significant
#' or not,the analysis of variances (anova) is conducted in order to
#' obtain f statistics with degree of freedom 1 and n-k-1. \code{tol}
#' is the tolerance quantity, usually 0.001, meaning, \code{f_local_search}
#' function will return the critical value of the f statistics correct
#' up to 3-decimal place. \code{is_echo} is the boolean variable,
#' if \code{is_echo} is FALSE, then \code{f_local_search} will return
#' numeric value only, and if \code{is_echo} is TRUE, then
#' \code{f_local_search} will return iteration table.
#'
#' @return a numeric value or an iteration table
#'
#' @note In order to check wether the critical value obtained from
#' \code{f_local_search} function is true or not, pf function can
#' be used
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
#' @seealso If one need the critical value of z for dynamic level of
#' significance between 0 and 1, then on can use \code{z_local_search}
#' function. Similarly, for t and chi square statistics.
#'
#' @importFrom stats pf runif
#'
#' @examples
#' # Given the level of significance is 0.025, with
#' # degree of freedom for sample 1 is 21 and
#' # degree of freedom for sample 2 is 27, then
#' # the critical value of f statistics is computed
#' # as follows.
#' f_local_search(
#' level_of_significance = 0.025, dfreedom1 = 21,
#' dfreedom2 = 27, tol = 0.001, max_iteration = 50,
#' ) #[1] 2.236959
#' pf(2.236959, df1= 21, df2= 27) #[1] 0.9750136
#' 1 - pf(2.236959, df1= 21, df2= 27) #[1] 0.02498644
#'
#' f_local_search(
#' level_of_significance = 0.025, dfreedom1 = 21,
#' dfreedom2 = 27, tol = 0.001, max_iteration = 50,
#' is_echo = TRUE
#' ) # return an iteration table
#' @export
f_local_search <- function(level_of_significance, dfreedom1, dfreedom2, tol, max_iteration, li = 0, ui = 50, is_echo = FALSE) {
    given_area <- 1 - level_of_significance # a.k.a. confidence level

    # data gathering
    li_vect <- c()
    ui_vect <- c()
    f_val_vect <- c()
    areavect <- c()
    tol_vect <- c()

    iter <- 0
    current_tol <- abs(li - ui)
    while (current_tol > tol && iter <= max_iteration) {
        tol_vect <- append(tol_vect, current_tol) # gather the tolerance quantity

        iter <- iter + 1
        f_val_estimate <- runif(n = 1, min = li, max = ui)
        calculated_area <- pf(q = f_val_estimate, df1 = dfreedom1, df2 = dfreedom2)

        # data gathering
        li_vect <- append(li_vect, li)
        ui_vect <- append(ui_vect, ui)
        f_val_vect <- append(f_val_vect, f_val_estimate)
        areavect <- append(areavect, calculated_area)

        # update the tabu by reducing the interval of
        # pseudo random number generator
        ifelse(
            calculated_area > given_area,
            ui <- f_val_estimate, li <- f_val_estimate
        )
        current_tol <- abs(li - ui)
    }

    f_interpolate <- (li + ui) / 2 # interpolate the f value

    # data gathering
    li_vect <- append(li_vect, li)
    ui_vect <- append(ui_vect, ui)
    f_val_vect <- append(f_val_vect, f_interpolate)
    areavect <- append(areavect, pf(q = f_interpolate, df1 = dfreedom1, df2 = dfreedom2))
    tol_vect <- append(tol_vect, abs(li - ui))

    # create iteration table as return value if is_echo is TRUE
    iteration_table <- data.frame(
        "iteration" = seq_len(length(f_val_vect)),
        "lower_bound" = li_vect,
        "upper_bound" = ui_vect,
        "f_estimate" = f_val_vect,
        "cumulative_probability" = areavect,
        "absolute_error" = tol_vect
    )

    # return phase. If is_echo is false, then return only f_value
    ifelse(is_echo,
           return(iteration_table), return(f_val_vect[length(f_val_vect)]))
}
