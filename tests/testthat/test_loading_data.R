library(cgmisc)
library(GenABEL)
data(cgmisc_data)

test_that("dataset can be loaded", {
  expect_equal(class(data)[1], 'gwaa.data')
  expect_equal(nids(data), 207)
  expect_equal(nsnps(data), 174375)
})
