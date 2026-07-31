if (requireNamespace("spelling", quietly = TRUE)) {
  not_cran <- Sys.getenv("NOT_CRAN")
  if (is.na(match(tolower(not_cran), c("1", "yes", "true")))) {
    cat("All Done!\n")
  } else {
    pkg_dir <- list.files("../00_pkg_src", full.names = TRUE)
    if (!length(pkg_dir)) {
      check_dir <- dirname(getwd())
      if (grepl("\\.Rcheck$", check_dir)) {
        source_dir <- sub("\\.Rcheck$", "", check_dir)
        if (file.exists(source_dir)) pkg_dir <- source_dir
      }
    }
    if (!length(pkg_dir) && identical(basename(getwd()), "tests")) {
      if (file.exists("../DESCRIPTION")) pkg_dir <- dirname(getwd())
    }
    if (!length(pkg_dir)) pkg_dir <- "."

    # spelling::spell_check_package(vignettes = TRUE) applies a single
    # package-wide dictionary (en-US, from DESCRIPTION) to every vignette
    # file, with no option to exclude individual files. The hand-translated
    # Japanese vignettes (vignettes/*-ja.Rmd, vignettes/articles/*-ja.Rmd)
    # are not English prose, so running them past an en-US hunspell
    # dictionary only produced a flood of false-positive word fragments
    # (see the git history of inst/WORDLIST). Check Rd/DESCRIPTION via the
    # normal package check, then separately check only the
    # English-language vignettes/README/NEWS.
    results_pkg <- spelling::spell_check_package(pkg_dir, vignettes = FALSE)

    md_files <- list.files(
      file.path(pkg_dir, "vignettes"),
      pattern = "\\.q?r?md$", ignore.case = TRUE,
      full.names = TRUE, recursive = TRUE
    )
    md_files <- c(md_files, list.files(
      pkg_dir, pattern = "(readme|news|changes|index).r?q?md",
      ignore.case = TRUE, full.names = TRUE
    ))
    md_files <- md_files[!grepl("-ja\\.rmd$", md_files, ignore.case = TRUE)]

    # Match spell_check_package()'s ignore list exactly (package name +
    # Authors@R names + hunspell's common-word list + WORDLIST), otherwise
    # this second pass would re-flag the author names in README.
    desc <- read.dcf(file.path(pkg_dir, "DESCRIPTION"))[1, ]
    author_words <- tryCatch(
      unlist(eval(parse(text = desc[["Authors@R"]])), recursive = TRUE, use.names = FALSE),
      error = function(e) character(0)
    )
    ignore_words <- unique(c(desc[["Package"]], author_words, hunspell::en_stats, spelling::get_wordlist(pkg_dir)))

    results_md <- spelling::spell_check_files(
      md_files,
      ignore = ignore_words,
      lang = "en-US"
    )

    if (nrow(results_pkg)) {
      cat("Potential spelling errors (Rd / DESCRIPTION):\n")
      print(results_pkg)
    }
    if (nrow(results_md)) {
      cat("Potential spelling errors (English vignettes / README / NEWS):\n")
      print(results_md)
    }
    if (nrow(results_pkg) || nrow(results_md)) {
      cat("If these are false positive, run `spelling::update_wordlist()`.")
    }
    cat("All Done!\n")
  }
}
