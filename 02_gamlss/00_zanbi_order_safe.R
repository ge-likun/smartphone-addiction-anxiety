ZANBI_ORDER_SAFE_ID <- paste0(
  "gamlss.dist::ZANBI order-safe zero-vector patch; gamlss.dist ",
  as.character(utils::packageVersion("gamlss.dist"))
)

zanbi_recycle <- function(...) {
  values <- list(...)
  target_length <- max(lengths(values))
  lapply(values, rep, length.out = target_length)
}

zanbi_validate_parameters <- function(mu, sigma, nu) {
  if (any(mu <= 0)) stop("mu must be greater than 0")
  if (any(sigma <= 0)) stop("sigma must be greater than 0")
  if (any(nu <= 0) || any(nu >= 1)) stop("nu must be between 0 and 1")
  invisible(TRUE)
}

dZANBI_order_safe <- function(x, mu = 1, sigma = 1, nu = 0.3, log = FALSE) {
  values <- zanbi_recycle(x = x, mu = mu, sigma = sigma, nu = nu)
  x <- values$x
  mu <- values$mu
  sigma <- values$sigma
  nu <- values$nu
  zanbi_validate_parameters(mu, sigma, nu)

  zero <- rep.int(0L, length(x))
  log_probability_zero_nbi <- gamlss.dist::dNBI(zero, mu = mu, sigma = sigma, log = TRUE)
  log_probability_nbi <- gamlss.dist::dNBI(x, mu = mu, sigma = sigma, log = TRUE)
  log_probability <- ifelse(
    x == 0,
    log(nu),
    log(1 - nu) + log_probability_nbi - log(1 - exp(log_probability_zero_nbi))
  )
  result <- if (isTRUE(log)) log_probability else exp(log_probability)
  ifelse(x < 0, 0, result)
}

pZANBI_order_safe <- function(q, mu = 1, sigma = 1, nu = 0.3,
                              lower.tail = TRUE, log.p = FALSE) {
  values <- zanbi_recycle(q = q, mu = mu, sigma = sigma, nu = nu)
  q <- values$q
  mu <- values$mu
  sigma <- values$sigma
  nu <- values$nu
  zanbi_validate_parameters(mu, sigma, nu)

  zero <- rep.int(0L, length(q))
  cdf_zero_nbi <- gamlss.dist::pNBI(zero, mu = mu, sigma = sigma)
  cdf_nbi <- gamlss.dist::pNBI(q, mu = mu, sigma = sigma)
  probability <- ifelse(
    q == 0,
    nu,
    nu + (1 - nu) * (cdf_nbi - cdf_zero_nbi) / (1 - cdf_zero_nbi)
  )
  if (!isTRUE(lower.tail)) probability <- 1 - probability
  if (isTRUE(log.p)) probability <- log(probability)
  ifelse(q < 0, 0, probability)
}

qZANBI_order_safe <- function(p, mu = 1, sigma = 1, nu = 0.3,
                              lower.tail = TRUE, log.p = FALSE) {
  if (isTRUE(log.p)) p <- exp(p)
  if (!isTRUE(lower.tail)) p <- 1 - p
  if (any(p < 0) || any(p > 1)) stop("p must be between 0 and 1")

  values <- zanbi_recycle(p = p, mu = mu, sigma = sigma, nu = nu)
  p <- values$p
  mu <- values$mu
  sigma <- values$sigma
  nu <- values$nu
  zanbi_validate_parameters(mu, sigma, nu)

  zero <- rep.int(0L, length(p))
  cdf_zero_nbi <- gamlss.dist::pNBI(zero, mu = mu, sigma = sigma)
  positive_probability <- (p - nu) / (1 - nu) - 1e-10
  nbi_probability <- cdf_zero_nbi * (1 - positive_probability) + positive_probability
  nbi_probability <- ifelse(nbi_probability > 0, nbi_probability, 0)
  gamlss.dist::qNBI(nbi_probability, mu = mu, sigma = sigma)
}

rZANBI_order_safe <- function(n, mu = 1, sigma = 1, nu = 0.3) {
  if (length(n) > 1L) n <- length(n)
  qZANBI_order_safe(stats::runif(n), mu = mu, sigma = sigma, nu = nu)
}

zanbi_order_safe_score_components <- function(y, mu, sigma) {
  values <- zanbi_recycle(y = y, mu = mu, sigma = sigma)
  y <- values$y
  mu <- values$mu
  sigma <- values$sigma
  zero <- rep.int(0L, length(y))
  probability_zero <- gamlss.dist::dNBI(zero, mu = mu, sigma = sigma)
  denominator <- 1 - probability_zero
  nbi <- gamlss.dist::NBI()
  list(
    mu = nbi$dldm(y, mu, sigma) +
      probability_zero * nbi$dldm(zero, mu, sigma) / denominator,
    sigma = nbi$dldd(y, mu, sigma) +
      probability_zero * nbi$dldd(zero, mu, sigma) / denominator
  )
}

zanbi_order_safe_dldm <- function(y, mu, sigma, nu) {
  score <- zanbi_order_safe_score_components(y, mu, sigma)$mu
  ifelse(y == 0, 0, score)
}

zanbi_order_safe_d2ldm2 <- function(y, mu, sigma, nu) {
  dldm <- zanbi_order_safe_dldm(y, mu, sigma, nu)
  -dldm * dldm
}

zanbi_order_safe_dldd <- function(y, mu, sigma, nu) {
  sigma <- ifelse(sigma < 1e-6, 1e-6, sigma)
  score <- zanbi_order_safe_score_components(y, mu, sigma)$sigma
  ifelse(y == 0, 0, score)
}

zanbi_order_safe_d2ldd2 <- function(y, mu, sigma, nu) {
  dldd <- zanbi_order_safe_dldd(y, mu, sigma, nu)
  d2ldd2 <- -dldd^2
  ifelse(d2ldd2 < -1e-10, d2ldd2, -1e-10)
}

zanbi_order_safe_d2ldmdd <- function(y, mu, sigma, nu) {
  sigma <- ifelse(sigma < 1e-6, 1e-6, sigma)
  scores <- zanbi_order_safe_score_components(y, mu, sigma)
  dldm <- ifelse(y == 0, 0, scores$mu)
  dldd <- ifelse(y == 0, 0, scores$sigma)
  -dldm * dldd
}

zanbi_order_safe_gdev <- function(y, mu, sigma, nu, ...) {
  -2 * dZANBI_order_safe(y, mu = mu, sigma = sigma, nu = nu, log = TRUE)
}

ZANBI_order_safe <- function(mu.link = "log", sigma.link = "log", nu.link = "logit") {
  family <- do.call(
    gamlss.dist::ZANBI,
    list(mu.link = mu.link, sigma.link = sigma.link, nu.link = nu.link)
  )
  family$dldm <- zanbi_order_safe_dldm
  family$d2ldm2 <- zanbi_order_safe_d2ldm2
  family$dldd <- zanbi_order_safe_dldd
  family$d2ldd2 <- zanbi_order_safe_d2ldd2
  family$d2ldmdd <- zanbi_order_safe_d2ldmdd
  family$G.dev.incr <- zanbi_order_safe_gdev
  family$rqres <- expression(
    rqres(
      pfun = "pZANBI_order_safe", type = "Discrete", ymin = 0,
      y = y, mu = mu, sigma = sigma, nu = nu
    )
  )
  attr(family, "zanbi_implementation") <- ZANBI_ORDER_SAFE_ID
  family
}

distribution_function_for <- function(family_id, prefix = c("d", "p", "q", "r")) {
  prefix <- match.arg(prefix)
  if (identical(family_id, "ZANBI")) {
    return(get(paste0(prefix, "ZANBI_order_safe"), mode = "function", inherits = TRUE))
  }
  get(paste0(prefix, family_id), envir = asNamespace("gamlss.dist"))
}

zanbi_order_safe_audit <- function(tolerance = 1e-10) {
  n <- 88L
  y <- rep(0:21, length.out = n)
  mu <- exp(seq(log(0.35), log(14), length.out = n))
  sigma <- seq(0.04, 1.25, length.out = n)
  nu <- seq(0.03, 0.88, length.out = n)
  probability <- seq(0.005, 0.995, length.out = n)
  permutation <- order((seq_len(n) * 37L) %% n)
  inverse_permutation <- order(permutation)

  maximum_difference <- function(first, second) {
    if (length(first) != length(second) || any(!is.finite(first)) || any(!is.finite(second))) {
      return(Inf)
    }
    max(abs(first - second))
  }
  scalar_evaluate <- function(fun, arguments) {
    vapply(seq_len(n), function(index) {
      do.call(fun, lapply(arguments, function(value) value[index]))
    }, numeric(1))
  }
  audit_function <- function(component, safe_fun, reference_fun, arguments) {
    vector_value <- do.call(safe_fun, arguments)
    scalar_value <- scalar_evaluate(reference_fun, arguments)
    reversed_value <- rev(do.call(safe_fun, lapply(arguments, rev)))
    permuted_value <- do.call(
      safe_fun,
      lapply(arguments, function(value) value[permutation])
    )[inverse_permutation]
    scalar_difference <- maximum_difference(vector_value, scalar_value)
    reverse_difference <- maximum_difference(vector_value, reversed_value)
    permutation_difference <- maximum_difference(vector_value, permuted_value)
    data.frame(
      component = component,
      scalar_max_abs_difference = scalar_difference,
      reverse_max_abs_difference = reverse_difference,
      permutation_max_abs_difference = permutation_difference,
      tolerance = tolerance,
      passed = all(c(scalar_difference, reverse_difference, permutation_difference) <= tolerance),
      stringsAsFactors = FALSE
    )
  }

  reference_family <- gamlss.dist::ZANBI()
  safe_family <- ZANBI_order_safe()
  distribution_tests <- list(
    audit_function(
      "density",
      dZANBI_order_safe,
      gamlss.dist::dZANBI,
      list(x = y, mu = mu, sigma = sigma, nu = nu)
    ),
    audit_function(
      "cdf",
      pZANBI_order_safe,
      gamlss.dist::pZANBI,
      list(q = y, mu = mu, sigma = sigma, nu = nu)
    ),
    audit_function(
      "quantile",
      qZANBI_order_safe,
      gamlss.dist::qZANBI,
      list(p = probability, mu = mu, sigma = sigma, nu = nu)
    )
  )
  derivative_arguments <- list(y = y, mu = mu, sigma = sigma, nu = nu)
  derivative_tests <- lapply(
    c("dldm", "d2ldm2", "dldd", "d2ldd2", "d2ldmdd", "G.dev.incr"),
    function(component) {
      audit_function(
        component,
        safe_family[[component]],
        reference_family[[component]],
        derivative_arguments
      )
    }
  )
  do.call(rbind, c(distribution_tests, derivative_tests))
}
