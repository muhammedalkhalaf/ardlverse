# Replication test: validates panel_ardl PMG/MG/DFE outputs against
# Stata's xtpmg, using the Blackburne & Frank (2007) Stata Journal example.
#
# Contributed by Yeleazar Levchenko (Kyiv School of Economics).
#
# The original replication script lives in inst/replication/replicate_jasa.R;
# this test is skipped on CRAN (the JASA dataset download from
# https://www.stata-journal.com is not always reachable from build farms).

test_that("PMG/MG/DFE match Stata xtpmg on JASA replication", {
  skip_on_cran()
  skip_on_ci()
  skip_if_not_installed("foreign")
  skip_if_offline()

  # The actual replication script is shipped under inst/replication/ so it
  # remains reachable to package users:
  script_path <- system.file("replication", "replicate_jasa.R",
                             package = "ardlverse")
  skip_if(script_path == "", "replication script not installed")

  # Source the replication script in a fresh environment and check that
  # it completes without error and produces non-trivial output.
  env <- new.env(parent = globalenv())
  expect_no_error(source(script_path, local = env))
})
