# Formula parsing shared by boot_ardl(), aardl(), mtnardl() and fbnardl().
# Only plain variable names are accepted: earlier versions passed the formula
# through all.vars(), which silently turned log(y) ~ log(x) into y ~ x.
.ardl_formula_vars <- function(formula, data) {
  if (!inherits(formula, "formula") || length(formula) != 3)
    stop("'formula' must be a two-sided formula such as y ~ x1 + x2", call. = FALSE)
  if (!is.data.frame(data)) stop("'data' must be a data frame", call. = FALSE)
  lhs <- formula[[2]]
  if (!is.name(lhs))
    stop("transformed terms such as log(y) are not supported in the formula; ",
         "create the transformed variable in 'data' and use its name",
         call. = FALSE)
  tt <- stats::terms(formula, data = data)
  labs <- attr(tt, "term.labels")
  rhs_vars <- all.vars(formula[[3]])
  if ("." %in% rhs_vars) rhs_vars <- labs
  if (length(labs) == 0)
    stop("the formula needs at least one regressor", call. = FALSE)
  if (!setequal(labs, rhs_vars) || any(!vapply(labs, function(l) {
    e <- str2lang(l)
    is.name(e)
  }, logical(1))))
    stop("transformed terms or interactions such as log(x), I(x^2) or x:z are ",
         "not supported in the formula; create the variables in 'data' and ",
         "use their names", call. = FALSE)
  y_var <- as.character(lhs)
  vars <- c(y_var, labs)
  miss <- setdiff(vars, names(data))
  if (length(miss))
    stop("variable(s) not found in 'data': ", paste(miss, collapse = ", "),
         call. = FALSE)
  for (v in vars) {
    if (!is.numeric(data[[v]]))
      stop("variable '", v, "' must be numeric", call. = FALSE)
    if (anyNA(data[[v]]))
      stop("variable '", v, "' has missing values; the time series must be ",
           "complete and without gaps", call. = FALSE)
  }
  list(y_var = y_var, x_vars = labs)
}
