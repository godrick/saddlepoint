test_that("testCGF_K3K4ops passes on representative base and wrapper cases", {
  cases <- make_testCGF_K3K4ops_cases()

  for (case_name in names(cases)) {
    res <- suppressMessages(do.call(testCGF_K3K4ops, cases[[case_name]]))
    expect_testCGF_K3K4ops_passes(case_name, res)
  }
})
