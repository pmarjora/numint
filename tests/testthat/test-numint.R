context("Basic set of tests for numint")
test_that("First test", {

  # Do a Normal test

  for (i in 1:1){
  ans <- num_int(dnorm, a = -1, b = 1, mean = 0, sd = 1,
                 N = 1e5, ncores = 1)
  ans
 ExpectedAns<- 1 - pnorm(-1)*2
  #plot(ans)

  # Is the answer correct
  expect_equal(ans$val,ExpectedAns,tolerance=0.01)
  }
})

test_that("Second test",{
  ans <- num_int(dnorm, a = -1, b = 1, mean = 0, sd = 1,
                 N = 1e5, ncores = 1)
  expect_s3_class(ans,"numint")
  expect_silent(plot(ans))
})
