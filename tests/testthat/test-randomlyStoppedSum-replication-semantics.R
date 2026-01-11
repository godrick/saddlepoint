test_that("RSS replication semantics: vector-valued Y vs i.i.d. replicates of scalar Y", {
  # This test illustrates the key ambiguity that motivates requiring block_size and/or iidReps:
  #
  # 1) "One replicate of a *vector*-valued Y"  (dimension = 2 here)
  #    -> coordinates share the SAME random count N, so cross-covariances can be nonzero.
  #
  # 2) "Two i.i.d. replicates of a *scalar* Y"
  #    -> replicates have independent counts N_1 and N_2, so cross-covariances are zero
  #       (block diagonal K2).

  count_cgf <- GeometricModelCGF(prob = adaptor(fixed_param = 0.3), iidReps = 1)

  # scalar Bernoulli summand (1D X)
  summand_cgf_scalar <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(fixed_param = 0.2),
    iidReps = 1
  )

  # vectorized Bernoulli summand: feeding a length-2 tvec means "2D X" (independent coords)
  summand_cgf_vector <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(fixed_param = 0.2),
    iidReps = "any"
  )

  t2 <- c(0.1, -0.2)
  theta <- 1  # unused because we used fixed_param adaptors

  # --- Case A: ONE 2D observation Y (shared N across coordinates) ---
  # We force "one replicate" using iidReps=1, and declare block_size=2 so tvec is interpreted as a single 2D vector.
  rss_vectorY <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf_vector,
    block_size = 2,
    iidReps = 1
  )

  K2_vec <- as.matrix(rss_vectorY$K2(t2, theta))
  expect_true(abs(K2_vec[1, 2]) > 1e-10)  # generally nonzero due to shared N

  # --- Case B: TWO i.i.d. scalar observations Y_1, Y_2 (independent N_1, N_2) ---
  # Here block_size=1, iidReps="any" means "split tvec into scalar blocks"; with length(tvec)=2 we get 2 blocks.
  rss_iid_scalarY <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_cgf_scalar,
    block_size = 1,
    iidReps = "any"
  )

  K2_iid <- as.matrix(rss_iid_scalarY$K2(t2, theta))
  expect_equal(K2_iid[1, 2], 0, tolerance = 1e-12)
  expect_equal(K2_iid[2, 1], 0, tolerance = 1e-12)
})



test_that("RSS replication semantics: multinomial summand example", {
  # This test mirrors the scalar example, but with a 3D multinomial summand.
  # It shows that:
  #  - If you build ONE 6D Y (block_size=6, iidReps=1), then cross-block terms can be nonzero
  #    because both 3D chunks share the same N.
  #  - If you build TWO i.i.d. 3D replicates of Y (block_size=3, iidReps='any'), cross-block terms are zero.

  set.seed(1)

  count_cgf <- PoissonModelCGF(lambda = adaptor(fixed_param = 2.0), iidReps = 1)

  summand_mn <- MultinomialModelCGF(
    n = adaptor(fixed_param = 1),
    prob_vec = adaptor(fixed_param = c(0.2, 0.3, 0.5)),
    iidReps = "any"
  )

  t6 <- rnorm(6) * 0.1
  theta <- 1  # unused because we used fixed_param adaptors

  # ONE replicate of a 6D Y (internally, summand interprets length-6 as 2 blocks of 3)
  rss_one6D <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_mn,
    block_size = 6,
    iidReps = 1
  )

  # TWO i.i.d. replicates of a 3D Y (outer replication)
  rss_two3D <- randomlyStoppedSumCGF(
    count_cgf = count_cgf,
    summand_cgf = summand_mn,
    block_size = 3,
    iidReps = "any"
  )

  K2_one6D <- as.matrix(rss_one6D$K2(t6, theta))
  K2_two3D <- as.matrix(rss_two3D$K2(t6, theta))

  # Cross-block indices: first 3 vs last 3
  cross_one6D <- K2_one6D[1:3, 4:6]
  cross_two3D <- K2_two3D[1:3, 4:6]

  expect_true(max(abs(cross_one6D)) > 1e-10)
  expect_equal(max(abs(cross_two3D)), 0, tolerance = 1e-12)
})



test_that("RSS constructor rejects ambiguous replication specs", {
  # The package uses vectorized summand CGFs heavily. For RSS, that can silently change model meaning.
  # So the constructor should reject calls where BOTH block_size and a numeric iidReps are missing.

  count_cgf <- PoissonModelCGF(lambda = adaptor(fixed_param = 2.0), iidReps = 1)
  summand_cgf <- BinomialModelCGF(
    n = adaptor(fixed_param = 1),
    p = adaptor(fixed_param = 0.3),
    iidReps = 1
  )

  ## expect error
  ## randomlyStoppedSumCGF(count_cgf, summand_cgf)

  # iidReps='any' without block_size is also ambiguous (cannot split tvec).
  expect_error(
    randomlyStoppedSumCGF(count_cgf, summand_cgf, iidReps = "any"),
    "block_size")
})
