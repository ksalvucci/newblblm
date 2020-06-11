m = 3
B = 1000
fit <- newblblm(mpg~wt, mtcars, m = m, B = B)

test_that("correct class", {
  expect_equal(class(fit), "newblblm")
})


test_that("correct subsets", {
  expect_equal(length(fit$estimates), m)
})


b <- NULL
for (i in fit$estimates) {
  num_boots <- length(i)
  b <- c(b, num_boots)
}

test_that("correct bootstraps", {
  expect_equal(b, rep(B,m))
})


lm1coef <- lm1(mpg ~ wt, mtcars, freq = NULL)$coef
lmcoef <- lm(mpg ~ wt, mtcars)$coefficients

test_that("lm works", {
  expect_equal(lm1coef, lmcoef)
})
