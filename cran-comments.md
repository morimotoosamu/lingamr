## Submission

This is a new submission, version 0.1.1, following version 0.1.0. It
carries over the previous round's fixes (software names in single quotes
in Title and Description; `\value` tags added to all exported functions;
`\examples` added to the exported `print.*` methods; unconditional
Suggests-package dependencies removed from examples; the "no visible
binding for global variable" NOTE fixed) and additionally:

* Fixed several correctness issues found during an internal review: a
  condition in the kernel-based independence measure that silently ignored
  soft prior knowledge, an unsupported `reg_method = "ridge"` combination in
  the bootstrap total-effect step, and a data-scale dependence in the
  default adaptive-LASSO regularization path. Test suite expanded
  accordingly.
* The kernel-based independence measure (`measure = "kernel"`) now uses an
  incomplete-Cholesky low-rank approximation for `n > 1000`, substantially
  reducing runtime and memory for large samples; behavior for `n <= 1000`
  is unchanged.
* Various smaller fixes and additions; see NEWS.md for the full list.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* Local: Windows 11, R 4.6.0

## Downstream dependencies

There are currently no downstream dependencies for this package.
