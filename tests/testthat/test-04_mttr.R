test_that("MTTR and its variance agree with two-state analytic results", {
  p <- 0.2
  q <- 0.3
  P <- matrix(c(1 - p, p, q, 1 - q), 2, byrow = TRUE)
  for (down_mass in c(1, 0.4)) {
    x <- MarkovChain(c(1 - down_mass, down_mass), P)
    expect_equal(mttr(x, 2), down_mass / q)
    expected <- down_mass^2 * (1 - q) * (p + q) / (p * q^3)
    expect_equal(var_mttr(x, 2), expected)
  }
})

test_that("MTTR wrappers preserve state selection and MTTF results", {
  labels <- c("failed_a", "working", "failed_b")
  P <- matrix(c(0.6, 0.3, 0.1, 0.2, 0.5, 0.3, 0.4, 0.1, 0.5),
              3, byrow = TRUE, dimnames = list(labels, labels))
  x <- MarkovChain(setNames(c(0.7, 0.2, 0.1), labels), P)
  for (down in list(c(1, 3), c(3, 1), labels[c(3, 1)])) {
    expect_equal(mttr(x, down), mttf(x, c(1, 3)))
    expect_equal(var_mttr(x, down), var_mttf(x, c(1, 3)))
  }
})

test_that("MTTR wrappers retain input validation", {
  x <- MarkovChain(c(0, 1), matrix(c(0.8, 0.2, 0.3, 0.7), 2,
                                byrow = TRUE))
  for (f in list(mttr, var_mttr)) {
    expect_error(f(x, c(2, 2)), "distinct")
    expect_error(f(x, integer()), "nonempty")
    expect_error(f(x, "unknown"), "unknown")
    expect_error(f(x, 1:2), "proper subset")
    invalid <- x
    invalid$a[1] <- -1
    expect_error(f(invalid, 2), "valid StochVec")
  }
})

test_that("MTTR allows transient down states but variance needs irreducibility", {
  x <- MarkovChain(c(0, 1), matrix(c(1, 0, 0.3, 0.7), 2, byrow = TRUE))
  expect_equal(mttr(x, 2), 1 / 0.3)
  expect_error(var_mttr(x, 2), "irreducible")
})
