# Try all tests
testthat::test_dir("tests")

# Alternative usage
testthat::test_file("./tests/test_nearest.R")
testthat::test_file("./tests/test_hole.R")
