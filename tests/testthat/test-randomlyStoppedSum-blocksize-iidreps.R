test_that("randomlyStoppedSumCGF: block_size and iidReps semantics (constructor equivalences)", {
  # This test checks that the three common ways of specifying replication are equivalent:
  #
  # (A) specify block_size only   -> iidReps inferred from length(tvec)  (i.e., iidReps='any')
  # (B) specify iidReps only      -> block_size inferred from length(tvec)/iidReps
  # (C) specify both explicitly

  count_cgf <- PoissonModelCGF(lambda = adaptor(indices = 1), iidReps = 1)
  summand_cgf <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(indices = 2),
    iidReps = 1
  )

  theta <- c(2.0, 0.3)  # lambda_N, p

  # 3 scalar observations => tvec length 3
  tvec <- c(0.10, -0.20, 0.05)

  # A) block_size only (iidReps inferred, i.e. "any")
  rss_A <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf,
    block_size = 1
    # iidReps omitted on purpose: should behave like iidReps='any'
  )

  # B) iidReps only (block_size inferred from length(tvec)/iidReps)
  rss_B <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf,
    iidReps = 3
  )

  # C) both
  rss_C <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf,
    block_size = 1,
    iidReps = 3
  )

  # Check K, K1, K2 match across constructors
  expect_equal(as.numeric(rss_A$K(tvec, theta)), as.numeric(rss_B$K(tvec, theta)), tolerance = 1e-12)
  expect_equal(as.numeric(rss_A$K(tvec, theta)), as.numeric(rss_C$K(tvec, theta)), tolerance = 1e-12)

  expect_equal(as.numeric(rss_A$K1(tvec, theta)), as.numeric(rss_B$K1(tvec, theta)), tolerance = 1e-12)
  expect_equal(as.numeric(rss_A$K1(tvec, theta)), as.numeric(rss_C$K1(tvec, theta)), tolerance = 1e-12)

  expect_equal(as.matrix(rss_A$K2(tvec, theta)), as.matrix(rss_B$K2(tvec, theta)), tolerance = 1e-12)
  expect_equal(as.matrix(rss_A$K2(tvec, theta)), as.matrix(rss_C$K2(tvec, theta)), tolerance = 1e-12)
})
