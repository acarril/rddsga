# rddsga
Subgroup analysis for regression discontinuity designs using inverse propensity score weighting

# Installation

The package can be installed directly from SSC:
```stata
ssc install rddsga
```

Alternatively, you can get the latest version hosted in this repository by first setting the location with `net`:
```stata
net from https://raw.githubusercontent.com/acarril/rddsga/master/
```
Then you can use the following commands to describe and install the package and its ancillary files:
```stata
// Describe
net describe rddsga
// Install ado-files and help files
net install rddsga
// Install ancillary files (datasets)
net get rddsga
```

If you have installed the program and want to update it, be sure to `net install` and `net get` using options `replace force`.

# Help

Refer to the included help file for details regarding usage.
For additional details regarding the methodology refer to the [project wiki](https://gitlab.com/acarril/rddsga/wikis/home).
