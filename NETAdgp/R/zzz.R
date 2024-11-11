# NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT Proprietary Software
# © 2024 NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT. All rights reserved.
#
# This software is licensed under the NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT Software License Agreement.
# You may not use, modify, or distribute this software except in compliance with the license.
# Please refer to the LICENSE file for the full terms of the license.


.onLoad <- function(libname, pkgname) {
  packageStartupMessage("Welcome to AutomationND. This software is licensed under the NEURAL EPILEPTOTRANSCRIPTOMIC ADVANCEMENT Proprietary License.")
  packageStartupMessage("Please see the LICENSE file for the full terms of the license.")
  required_packages <- c("dplyr", "clustermole", "Seurat", "STRINGdb", "tidyr", 
                        "limma", "gprofiler2", 
                         "edgeR")
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      install.packages(pkg)
    }
  }
}
