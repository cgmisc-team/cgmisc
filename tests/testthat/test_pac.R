library(cgmisc)
library(GenABEL)

test_that("pop allele counts work", {
  data.qc1 <- readRDS('data.qc1.rds')
  pop <- readRDS('pop.rds')
  pac <- pop.allele.counts(data = data.qc1[ ,data.qc1@gtdata@chromosome == 2], pops = pop, progress=F)
  expect_equal_to_reference(pac, 'pac.rds')
})