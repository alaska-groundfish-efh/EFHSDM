# Scripts use to run pkg development

# Documentation
devtools::document()


# Dependencies - need to add everything that gets called w library() or require() here.
use_package(package, type = "Imports", min_version = NULL)
