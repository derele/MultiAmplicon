[![Travis-CI Build Status](https://travis-ci.org/derele/MultiAmplicon.svg?branch=master)](https://travis-ci.org/derele/MultiAmplicon) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/derele/MultiAmplicon?branch=master&svg=true)](https://ci.appveyor.com/project/derele/MultiAmplicon) [![Coverage status](https://codecov.io/gh/derele/MultiAmplicon/branch/master/graph/badge.svg)](https://codecov.io/github/derele/MultiAmplicon?branch=master)

# The MultiAmplicon R-package

The MultiAmplicon-package allows the matching and removal of primer
sequences from sequencing reads. Amplicons, here defined as the
sequences amplified by particular combination of forward and reverse
primer, can then be further processed. The package guiedes downstream
analysis of multiple amplicons in a bioconductor workflow wrapped
around the package
[dada2](https://benjjneb.github.io/dada2/index.html) and helps to
export to [phyloseq](https://joey711.github.io/phyloseq/index.html)
for further analysis.


## Installation
```S
require(devtools)
devtools::install_github("derele/MultiAmplicon")
```

## Documentation

Documentation, including the
[reference manual](https://derele.github.io/MultiAmplicon/reference/index.html)
is available
[online](https://derele.github.io/MultiAmplicon/index.html).

Tutorials cover:

- [A a typical workflow quickly executed on data available within the package](https://derele.github.io/MultiAmplicon/articles/MultiAmplicon-small-example.html)

- [A a typical workflow quickly executed on a real world dataset](https://derele.github.io/MultiAmplicon/articles/MultiAmplicon-real-world-example.html)


