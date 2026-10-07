library(testthat)
library(MMEFM)

# The bootstrap has no tests yet; run the suite once test files are added.
if (length(list.files("testthat", pattern = "^test.*[.]R$")) > 0L) {
    test_check("MMEFM")
}
