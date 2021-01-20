
#' Runs the analytic workflow for the Neuro project
#'
#' @keywords 4CE
#' @export

runAnalysis <- function() {

    ## make sure this instance has the latest version of the quality control and data wrangling code available
    devtools::install_github("https://github.com/covidclinical/Phase2.1DataRPackage", subdir="FourCePhase2.1Data", upgrade=FALSE)

    ## get the site identifier assocaited with the files stored in the /4ceData/Input directory that 
    ## is mounted to the container
    currSiteId = FourCePhase2.1Data::getSiteId()

    ## run the quality control
    FourCePhase2.1Data::runQC(currSiteId)

    ## DO NOT CHANGE ANYTHING ABOVE THIS LINE

    ## To Do: implement analytic workflow, saving results to a site-specific 
    ## file to be sent to the coordinating site later via submitAnalysis()

    ## Save results to appropriately named files for submitAnalysis(), e.g.:
    #write.csv(
    #    matrix(rnorm(100), ncol=5), 
    #    file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_ResultTable.csv"))
    #)

    #write.table(
    #    matrix(rnorm(12), ncol=3), 
    #    file=file.path(getProjectOutputDirectory(), paste0(currSiteId, "_ModelParameters.txt"))
    #)
    
}

