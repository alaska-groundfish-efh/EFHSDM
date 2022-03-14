# Scripts use to run pkg development


# Documentation -----------------------------------------------------------
# Package-level documentation
use_package_doc(open = rlang::is_interactive())

# Function documentation
devtools::document()

# Vignettes
usethis::use_vignette("EFHSDM")
devtools::build_vignettes()


# Dependencies ------------------------------------------------------------
# ***Need to add everything that gets called w library() or require() here.
use_package(package, type = "Imports", min_version = NULL)


# Datasets ----------------------------------------------------------------

use_data() #FILL IN with example datasets that need to be stored with the pkg
