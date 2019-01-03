sessionInfo()

library(devtools)
install_version("sp", version = "1.2-3", repos = "http://cran.us.r-project.org")
install_version("rgdal", version = "1.1-10", repos = "http://cran.us.r-project.org")

library(sp)
library(rgdal)
sessionInfo()