#' Finding the Value of z Statistics for Dynamic Significance Level
#'
#' \code{z_local_search} is used to find the critical value of z
#' statistics at level of significance between 0 and 1
#'
#' @usage z_local_search(level_of_significance, tol, max_iteration,
#' li = -7, ui = 7, is_echo = FALSE)
#'
#' @param level_of_significance a value between 0 and 1
#' @param tol the tolerance quantity, usually 0.001
#' @param max_iteration the maximum iteration as termination criteria of searching
#' @param li the lower bound of the searching
#' @param ui the upper bound of the search space
#' @param is_echo a boolean value which is by default, it's FALSE.
#'
#' @details Function z_local_search is used to find the critical
#' value of z statistics at level of significance between 0 and 1.
#' The critical value is used in the construction of interval to
#' estimate a certain parameter, and it's used in hypothesis testing
#' as well.
#'
#' Level of significance is correspond to the confidence limits, if
#' confidence limits is known, then the level of significance can be
#' obtained by subtracting the value of confidence limits from 1.
#' Tolerance quantity is used to find this critical point statistics
#' correct up to i-th decimal places, usually it's set to be 0.001.
#' The tolerance quantity is used as terminaton criteria too.
#' Maximum iteration max_itear is used to the termination criteria
#' too. The li and ui are, respectively, the lower bound and upper
#' bound of the searching. For z-statistics, by default, the value
#' of li and ui are, respectively, -7 and 7. The value of is_echo
#' is set ot be FALSE by default, since most of people only want to
#' obtained this critical value. If is_echo is TRUE, then the
#' function will return an iteration table as data frame.
#'
#' @return a number (exclusive) or data frame.
#'
#' @note In order to check wether the return value is indeed, the critical value,
#' use pnorm function
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
#' @seealso If one needs critical value of t, then use t_local_search
#' function. Similarly for chisq_local_search, and f_local_search
#'
#' @importFrom stats pnorm runif
#'
#' @examples
#' # Given the level of significance is 0.025,
#' # then the critical value for z at this level of
#' # significance is computed as follows.
#' z_local_search(level_of_significance = 0.025, tol = 0.001,
#' max_iteration = 50
#' ) #[1] 1.960062
#' pnorm(1.960062) #[1] 0.9750057
#' 1-pnorm(1.960062) #[1] 0.02599427 ~= 0.025
#'
#' z_local_search(level_of_significance = 0.025, tol = 0.001,
#' max_iteration = 50, is_echo = TRUE) # return a data frame
#' @export
z_local_search <- function(level_of_significance, tol, max_iteration, li = -7, ui = 7, is_echo = FALSE) {
    given_area <- 1 - level_of_significance # a.k.a. confidence level

    # data gathering
    li_vect <- c()
    ui_vect <- c()
    zvect <- c()
    areavect <- c()
    tol_vect <- c()

    iter <- 0
    current_tol <- abs(li - ui)
    while (current_tol > tol && iter <= max_iteration) {
        tol_vect <- append(tol_vect, current_tol) # gather the tolerance quantity

        iter <- iter + 1
        zval_estimate <- runif(n = 1, min = li, max = ui)
        calculated_area <- pnorm(zval_estimate)

        # data gathering
        li_vect <- append(li_vect, li)
        ui_vect <- append(ui_vect, ui)
        zvect <- append(zvect, zval_estimate)
        areavect <- append(areavect, calculated_area)

        # update the tabu by reducing the interval of
        # pseudo random number generator
        ifelse(calculated_area > given_area,
               ui <- zval_estimate, li <- zval_estimate)

        current_tol <- abs(li - ui)
    }

    z_interpolate <- (li + ui)/2 # interpolate the z by using bisection

    # data gathering
    li_vect <- append(li_vect, li)
    ui_vect <- append(ui_vect, ui)
    zvect <- append(zvect, z_interpolate)
    areavect <- append(areavect, pnorm(z_interpolate))
    tol_vect <- append(tol_vect, abs(li - ui))

    # create iteration table as return value if is_echo == TRUE
    iteration_table <- data.frame(
        "iteration" = seq_len(length(zvect)),
        "lower_bound" = li_vect,
        "upper_bound" = ui_vect,
        "z_estimate" = zvect,
        "cumulative_probability" = areavect,
        "absolute_error" = tol_vect
    )

    # return phase. If is_echo is false, then return only z_value
    ifelse(is_echo,
           return(iteration_table), return(zvect[length(zvect)]))
}
