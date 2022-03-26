# Contributing

## Welcome
Thank you for your interest in contributing to `EFHSDM`! We plan to maintain this package as a resource for stakeholders, modelers, and researchers interested in species distribution modeling and the description of Essential Fish Habitat.

## Table of Contents

[Short Links](#short-links)

[How Can I Contribute?](#how-can-i-contribute)

[Style Guides](#style-guides)


## Short Links
* View other EFH work products or request data products [here](https://github.com/alaska-groundfish-efh)
* Report bugs [here](https://github.com/alaska-groundfish-efh/EFHSDM/issues)
* Contact maintainers: Margaret at margaret (dot) siple (at) noaa.gov and Jeremy Harris, jeremy (dot) harris (at) noaa.gov.

## How Can I Contribute?
 * [Reporting Bugs](#reporting-bugs)
 * [Suggesting Enhancements](#suggesting-enhancements)
 * [Pull Requests](#pull-requests)
 
## Style Guides
* [Write an issue on GitHub](#issue-style-guides)
* [Code style guides](#code-style-guides)

## Reporting Bugs 
Bugs are tracked as [GitHub issues](https://guides.github.com/features/issues/). If you find a bug, please make your report as detailed as possible.

**Note:** If you find a [**Closed** issue](https://github.com/alaska-groundfish-efh/EFHSDM/issues?q=is%3Aissue+is%3Aclosed) that looks like it is the same thing that you're experiencing, open a new issue and include a link to the original issue in the body of your new one (you can tag issues by a hashtag (#) followed by the issue number).

## Suggesting Enhancements 
:pencil: We are always interested in how to make this package more useful and accessible. If you have an enhancement to suggest, please submit an issue and tag it as a feature request. 

### How do I submit a good enhancement suggestion?
Enhancement suggestions are tracked as [GitHub issues](https://guides.github.com/features/issues/). Create an issue on the `EFHSDM` repository and include the following:

* **Use a clear and descriptive title** for the issue to identify the suggestion.
* **Provide a step-by-step description of the suggested enhancement** in as many details as possible.
* **Provide specific examples to demonstrate the steps**. Include copy/pasteable snippets which you use in those examples, as [Markdown code blocks](https://help.github.com/articles/markdown-basics/#multiple-lines).
* **Describe the current behavior** and **explain which behavior you expected to see instead** and why.
* **Include screenshots and animated GIFs** which help you demonstrate the steps or point out the part of Atom which the suggestion is related to. You can use [this tool](https://www.cockos.com/licecap/) to record GIFs on macOS and Windows, and [this tool](https://github.com/colinkeenan/silentcast) or [this tool](https://github.com/GNOME/byzanz) on Linux.
* **Explain why this enhancement would be useful** to most `EFHSDM` users and stakeholders.
* **Specify the name and version of the OS you're using.**


## Pull requests
If you have a code suggestion in hand, please submit a pull request. 

## Write an issue on GitHub
### Git commit messages
(many of these are lifted from the Atom [style guide](https://github.com/atom/atom/blob/master/CONTRIBUTING.md) for commit messages)

* Use the present tense ("Add feature" not "Added feature")
* Use the imperative mood ("Move cursor to..." not "Moves cursor to...")
* Limit the first line to 72 characters or less
* Reference issues and pull requests liberally after the first line
* When only changing documentation, include `[ci skip]` in the commit title
* Consider starting the commit message with an applicable emoji:
    * :art: `:art:` when improving the format/structure of the code
    * :racehorse: `:racehorse:` when improving performance
    * :memo: `:memo:` when writing docs
    * :penguin: `:penguin:` when fixing something on Linux
    * :apple: `:apple:` when fixing something on macOS
    * :checkered_flag: `:checkered_flag:` when fixing something on Windows
    * :bug: `:bug:` when fixing a bug
    * :fire: `:fire:` when removing code or files
    * :white_check_mark: `:white_check_mark:` when adding tests
    * :arrow_up: `:arrow_up:` when upgrading dependencies
    * :arrow_down: `:arrow_down:` when downgrading dependencies


### Coding style
* Guides
  * Follow the [tidyverse style guide](https://style.tidyverse.org)
  * References on R package development from Hadley's book [R Packages](http://r-pkgs.had.co.nz)
  * [Advanced R](http://adv-r.had.co.nz)
  * [CRAN page](http://cran.r-project.org/doc/manuals/r-release/R-exts.html)
  * [Semantic versioning](https://semver.org/)
* Code rules
  * Use `vapply` or `lapply` instead of `sapply` because `sapply` is not consistent in the what class of object it returns (this may not yet be implemented in the package)
  * Use `[` rather than `subset` because `subset` will lead to unassigned objects, i.e., R will not know where the column vector name in the `subset` call comes from
  * Use `<-` rather than `=` to assign objects
  * Use `seq_len(NROW(x))` or `seq_along(x)` rather than `1:NROW(x)` or `1:length(x)` (not yet implemented in the package)
  * [Importing functions from other packages in ss3sim functions](Importing-functions-from-other-packages)

### Importing functions from other packages
(these guidelines are lifted and adapted from the ss3sim [developer wiki](https://github.com/ss3sim/ss3sim/wiki/developers), which contains a lot of useful resources for developing code)

What if you need to use a function from another package? Normally you might attach all the functions from a package with a call to `library()`. This is a bad practice in the context of writing an R package. We shouldn't be attaching other packages in the user's environment. In the context of writing a package we will instead import that function for our use. The package structure ensures that the user has the package installed.

For more details, see http://r-pkgs.had.co.nz/namespace.html

Say you want to use the function `mle2` from the bbmle package.

Here are the steps to use the function within your package:

* Add the package to the list of `Imports` in the `DESCRIPTION` file. If it's possible that functionality has changed in the past for the function then specify a `>=` version number for the package. You can find the version you have installed in many ways. One way is to type `help(package = "bbmle")` within R. You can also look at your `sessionInfo()` in R.
* In the function Rogxyen documentation include `@importFrom bbmle mle2`. If there are multiple functions to import than list them as `@importFrom packagename function1 function2` etc.
* Call the function as you would normally, e.g., `mle2(x, ...)`

The above steps outline one of two preferred notations. Within `EFHSDM` we also use the `::` notation to be explicit about what package a function comes from. E.g., `raster::raster()` clearly tells readers of the code that the `raster()` function comes from the `raster` package without readers having to navigate to the top of the code to read the roxygen documentation. This notation also requires that packages be listed in the `Imports` section of the `DESCRIPTION` file. We support either notation, but we have largely switched to this one as time goes on.

Occasionally your code might call so many functions from another package that you just want all the functions available. In that case use `@import bbmle` (for example) at the top. We have done this in a few places in EFHSDM. But it is best practices to import specific functions that we need because it saves time, which users appreciate. It also makes the code more explicit.


## Code of conduct
Contributors to `EFHSDM` follow the [Contributor Covenant](https://www.contributor-covenant.org/version/2/1/code_of_conduct/) code of conduct. For answers to common questions about this code of conduct, see https://www.contributor-covenant.org/faq.
