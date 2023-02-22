#' Finding the Value of t Statistics for Dynamic Significance Level
#'
#' \code{t_local_search} is used to find the critical value of
#' t statistics at level of significance between 0 and 1
#'
#' @usage t_local_search(level_of_significance, dfreedom, tol,
#' max_iteration, li, ui, is_echo)
#'
#' @param level_of_significance a value between 0 and 1
#' @param dfreedom the degree of freedom, it's positive integer
#' @param tol tolerance quantity
#' @param max_iteration maximum iteration of the searching
#' @param li lower bound of the search space
#' @param ui upper bound of the search space
#' @param is_echo a boolean value, by default, it's FALSE
#'
#' @details \code{t_local_search} function is used to find the
#' critical value of t at level of significance between 0 and 1.
#' It's used in the construction of 95\% confidence interval of
#' the mean parameter, and it's used in the hypothesis testing as
#' wells. This function contain 7 parameters, 4 parameter which need
#' to be filled, and 3 parameters which by default already filled.
#' These parameter are \code{level_of_significance}, \code{dfreedom},
#' \code{tol}, \code{max_iteration}, \code{li}, \code{ui}, and is_echo.
#'
#' \code{level_of_significance} is a value between 0 and 1,
#' dfreedom is the degree of freedom for given system. For example,
#' the degree of freedom for multilinear regression system is n-k-1
#' where k is the number of parameter. The \code{tol} parameter
#' is the tolerance quantity, usually it's set to be 0.001. If the
#' tolerance quantity is 0.001, then \code{t_local_search} function
#' will return the critical value of t statistics correct for
#' 3 decimal places. The \code{max_iteration} parameter is the
#' desired number iteration of the search. The tolerance quantity
#' \code{tol} and the maximum iteration \code{max_iteration}
#' are used as the termination criteria of the loop invariant.
#' The \code{li} and \code{ui} parameters are, respectively, the
#' lower bound and the upper bound of the search space. And
#' \code{is_echo} is parameter to determine either the
#' \code{t_local_search} function will return a numeric value or iteration table.
#'
#' @return Return a numeric value or iteration table
#'
#' @note In order to check either the return value of
#' \code{t_local_search} function is indeed, the critical value for
#' given level of significance or not, then the function \code{pt}
#' can be used
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
#' @seealso If one need the critical value of z statistics, then use
#' \code{z_local_search} function. Similarly, for chi square and f statistics.
#'
#' @importFrom stats pt runif
#'
#' @examples
#' # Given the level of significance is 0.025, and
#' # degree of freedom is n-1 = 29, then the
#' # critical value of t statistics is computed as follows.
#' t_local_search(
#' level_of_significance = 0.025, dfreedom = 29,
#' tol = 0.001, max_iteration = 50
#' ) #[1] 2.045238
#' pt(2.045238, df = 29) #[1] 0.950004
#' 1 - pt(2.045238, df = 29) #[1] 0.02499956 ~= 0.025
#' @export
t_local_search <- function(level_of_significance, dfreedom, tol, max_iteration, li = -7, ui = 7, is_echo = FALSE) {
    given_area <- 1 - level_of_significance # a.k.a confidence level

    # data gathering
    li_vect <- c()
    ui_vect <- c()
    t_val_vect <- c()
    areavect <- c()
    tol_vect <- c()

    iter <- 0
    current_tol <- abs(li - ui)
    while (current_tol > tol && iter <= max_iteration) {
        tol_vect <- append(tol_vect, current_tol) # gather the tolerance quantity

        iter <- iter + 1
        t_val_estimate <- runif(n = 1, min = li, max = ui)
        calculated_area <- pt(q = t_val_estimate, df = dfreedom)

        # data gathering
        li_vect <- append(li_vect, li)
        ui_vect <- append(ui_vect, ui)
        t_val_vect <- append(t_val_vect, t_val_estimate)
        areavect <- append(areavect, calculated_area)

        # update the tabu by reducing the interval of
        # pseudo random number generator
        ifelse(
            calculated_area > given_area,
            ui <- t_val_estimate, li <- t_val_estimate
        )
        current_tol <- abs(li - ui)
    }

    t_interpolate <- (li + ui) / 2 # interpolate the t value

    # data gathering
    li_vect <- append(li_vect, li)
    ui_vect <- append(ui_vect, ui)
    t_val_vect <- append(t_val_vect, t_interpolate)
    areavect <- append(areavect, pt(q = t_interpolate, df = dfreedom))
    tol_vect <- append(tol_vect, abs(li - ui))

    # create iteration table as return value if is_echo == TRUE
    iteration_table <- data.frame(
        "iteration" = seq_len(length(t_val_vect)),
        "lower_bound" = li_vect,
        "upper_bound" = ui_vect,
        "t_estimate" = t_val_vect,
        "cumulative_probability" = areavect,
        "absolute_error" = tol_vect
    )

    # return phase. If is_echo is false, then return only t_value
    ifelse(is_echo,
           return(iteration_table), return(t_val_vect[length(t_val_vect)]))
}
