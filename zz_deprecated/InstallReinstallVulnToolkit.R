# Uninstall and reinstall latest branch from GitHub

# 1. If package is loaded and in the memory, forget
if ("VulnToolkit" %in% (.packages())){
  detach("package:VulnToolkit", unload=TRUE) 
}

# 2. If remotes is not already installed, install it
if (! ("remotes" %in% installed.packages())) {
  install.packages("remotes")
}

# 3. Install package from developer branch of GitHub
devtools::install_github("https://github.com/Smithsonian/VulnToolkit")

# 4. Load version into memory
library(VulnToolkit)

