verbose = TRUE

if(verbose){
    startingTime = Sys.time()
    currentTime = startingTime
    cat("Loading libraries and setting up environments \n")
}

## Load the libraries
suppressMessages({
    library(faosws)
    library(faoswsAupus)
    library(faoswsUtil)
    library(data.table)
    library(igraph)
})

## Set up testing environments
if(Sys.getenv("USER") == "mk"){
    GetTestEnvironment(
        baseUrl = "https://hqlprswsas1.hq.un.fao.org:8181/sws",
        token = "1b5cb7c8-6acb-41df-bd73-fa3d42d01be6"
        )
}

## Setting up the fbs element codes
FBSelements =
    c("Value_measuredElementFS_51", "Value_measuredElementFS_61",
      "Value_measuredElementFS_91", "Value_measuredElementFS_101",
      "Value_measuredElementFS_111", "Value_measuredElementFS_121",
      "Value_measuredElementFS_141", "Value_measuredElementFS_151")


if(verbose){
    endTime = Sys.time()
    timeUsed = endTime - currentTime
    cat("\t Time used:",timeUsed, attr(timeUsed, "units") , "\n")
    currentTime = endTime
}

## Functions to save the data back
SaveAupusData = function(aupusData, nodes){
    updatedAupus = nodes[, colnames(aupusData), with = FALSE]
    updatedAupus[, timePointYearsSP := as.character(timePointYearsSP)]
    setnames(updatedAupus, "timePointYearsSP", "timePointYears")
    SaveData(domain = "faostat_one", dataset = "FS1_SUA",
             data = updatedAupus, normalize = FALSE)
}

SaveInputFromProcessingData = function(inputData, edges){
    updatedInput = edges[, colnames(inputData), with = FALSE]
    updatedInput[, timePointYearsSP := as.character(timePointYearsSP)]
    setnames(updatedInput, old = c("Value_input", "flagFaostat_input"),
             new = c("Value", "flagFaostat"))
    SaveData(domain = "faostat_one", dataset = "input_from_proc_fs",
             data = updatedInput, normalize = TRUE)
}

SavePopulationData = function(populationData){
    updatedPopulation = copy(populationData)
    updatedPopulation[, measuredItemFS := "1"]
    updatedPopulation[, timePointYearsSP := as.character(timePointYearsSP)]
    setnames(updatedPopulation, "timePointYearsSP", "timePointYears")
    setnames(updatedPopulation, old = colnames(updatedPopulation),
             new = gsub("population", "measuredElementFS",
                 colnames(updatedPopulation)))
    setcolorder(updatedPopulation, c("geographicAreaFS", "measuredItemFS",
                                     "timePointYears",
                                     "Value_measuredElementFS_11",
                                     "flagFaostat_measuredElementFS_11",
                                     "Value_measuredElementFS_21",
                                     "flagFaostat_measuredElementFS_21"))
    SaveData(domain = "faostat_one", dataset = "FS1_SUA",
             data = updatedPopulation, normalize = FALSE)
}

## Run the module by area
areaCodes = swsContext.datasets[[1]]@dimensions$geographicAreaFS@keys
## areaCodes = c("100", "231")
for(areas in areaCodes){
    if(verbose){
        cat("Running Aupus module for country", areas,
            "\n----------------------------------------------------------\n")
        cat("Setting parameters and getting datasets\n")
    }
    
    ## Get the parameter
    param = getAupusParameter(areaCode = areas, assignGlobal = FALSE)

    ## Get the data sets
    getAupusDataset(aupusParam = param)

    ## NOTE (Michael): This is a hack to fill in the missing columns
    missingColumns =
        c(paste0("Value_measuredElementFS_", c(541, 546)),
          paste0("flagFaostat_measuredElementFS_", c(541, 546)))
    aupusData[, `:=`(c(missingColumns),
                     list(as.numeric(NA), as.numeric(NA),
                          as.character(NA), as.character(NA)))]

    ## Run the whole aupus module and save the data back
    aupusModule =
        try(
            {
                if(verbose){
                    endTime = Sys.time()
                    timeUsed = endTime - currentTime
                    cat("\t Time used:",timeUsed, attr(timeUsed, "units") , "\n")
                    currentTime = endTime
                    cat("Running Aupus module\n")
                }

                ## Construct the aupus network representation
                aupusNetwork =
                    suaToNetworkRepresentation(extractionRateData =
                                                   extractionRateData,
                                               shareData = shareData,
                                               inputData = inputData,
                                               ratioData = ratioData,
                                               balanceElementData =
                                                   balanceElementData,
                                               itemInfoData = itemInfoData,
                                               populationData = populationData,
                                               aupusParam = param)

                ## Run the aupus to update the data
                updatedAupusNetwork =
                    with(aupusNetwork,
                         Aupus(nodes = nodes, edges = edges,
                               from = param$keyNames$itemParentName,
                               to = param$keyNames$itemChildName,
                               aupusParam = param))
                
                if(verbose){
                    endTime = Sys.time()
                    timeUsed = endTime - currentTime
                    cat("\t Time used:",timeUsed, attr(timeUsed, "units") , "\n")
                    currentTime = endTime
                    cat("Perform FBS standardization\n")
                }

                ## Construct the network for standardization
                standardizationGraph = 
                    with(updatedAupusNetwork,
                         constructStandardizationGraph(nodes = nodes,
                                                       edges = edges,
                                                       standardizeElement =
                                                           FBSelements,
                                                       from = param$keyNames$itemChildName,
                                                       to = param$keyNames$itemParentName,
                                                       aupusParam = param))

                ## Standardize the graph to get the FBS
                fbs =
                    fbsStandardization(graph = standardizationGraph,
                                       standardizeElement = FBSelements,
                                       aupusParam = param,
                                       plot = FALSE)
                
                if(verbose){
                    endTime = Sys.time()
                    timeUsed = endTime - currentTime
                    cat("\t Time used:",timeUsed, attr(timeUsed, "units") , "\n")
                    currentTime = endTime
                    cat("Save data back\n")
                }                
                ## Save the data back
                SaveAupusData(aupusData = aupusData,
                              nodes = updatedAupusNetwork$nodes)
                SaveInputFromProcessingData(inputData = inputData,
                                            edges = updatedAupusNetwork$edges)
                SavePopulationData(populationData = populationData)
                if(verbose){
                    endTime = Sys.time()
                    timeUsed = endTime - currentTime
                    cat("\t Time used:",timeUsed, attr(timeUsed, "units") , "\n")
                    currentTime = endTime
                }                
            }
        )
}

endTime = Sys.time()
totalTime = endTime - startingTime

## Return the results
if(inherits(aupusModule, "try-error")){
    cat("Aupus Module Failed after", totalTime, attr(totalTime, "units"), "\n")
} else {
    cat("Aupus Module Executed Successfully after", totalTime,
        attr(totalTime, "units"), "\n")
}
