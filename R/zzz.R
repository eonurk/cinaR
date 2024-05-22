#' @importFrom utils install.packages installed.packages
#' @importFrom BiocManager install
NULL

# Function to check and install missing Bioconductor dependencies
check_and_install_bioc_dependencies <- function() {
    missing_pkgs <- c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene", 
                      "TxDb.Hsapiens.UCSC.hg19.knownGene", 
                      "TxDb.Mmusculus.UCSC.mm10.knownGene")
    installed <- utils::installed.packages()[, "Package"]
    missing <- missing_pkgs[!(missing_pkgs %in% installed)]
    
    if (length(missing)) {
        message("Installing missing Bioconductor packages: ", paste(missing, collapse = ", "))
        if (!requireNamespace("BiocManager", quietly = TRUE)) {
            utils::install.packages("BiocManager")
        }
        BiocManager::install(missing)
    } else {
        message("All required Bioconductor packages are already installed.")
    }
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage("To ensure all required Bioconductor packages are installed, run: check_and_install_bioc_dependencies()")
}
