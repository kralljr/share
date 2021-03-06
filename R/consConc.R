#' Pollution concentrations over time for 85 northeast US monitors 
#' 
#' \itemize{
#' 	\item Date. Date in format %Y-%m-%d 
#'     	\item PM25_SPEC. Total PM2.5 mass measured at monitor (mug/m^3)
#'     	\item arsenic. Total PM2.5 measured at monitor (mug/m^3)
#'     	\item aluminum. Aluminum PM2.5 measured at monitor (mug/m^3)
#'     	\item bromine. Bromine PM2.5 measured at monitor (mug/m^3)
#'     	\item calcium. Calcium PM2.5 measured at monitor (mug/m^3)
#'     	\item copper. Copper PM2.5 measured at monitor (mug/m^3)
#'     	\item chlorine. Chlorine PM2.5 measured at monitor (mug/m^3)
#'     	\item iron. Iron PM2.5 measured at monitor (mug/m^3)
#'     	\item lead. Lead PM2.5 measured at monitor (mug/m^3)
#'     	\item manganese. Manganese PM2.5 measured at monitor (mug/m^3)
#'     	\item nickel. Nickel PM2.5 measured at monitor (mug/m^3)
#'     	\item phosphorus. Phosphorus PM2.5 measured at monitor (mug/m^3)
#'     	\item selenium. Selenium PM2.5 measured at monitor (mug/m^3)
#'     	\item titanium. Titanium PM2.5 measured at monitor (mug/m^3)
#'     	\item vanadium. Vanadium PM2.5 measured at monitor (mug/m^3)
#'     	\item silicon. Silicon PM2.5 measured at monitor (mug/m^3)
#'     	\item zinc. Zinc PM2.5 measured at monitor (mug/m^3)
#'     	\item strontium. Strontium PM2.5 measured at monitor (mug/m^3)
#'     	\item potassium. Potassium PM2.5 measured at monitor (mug/m^3)
#'     	\item ammonium_ion. Ammonium ion PM2.5 measured at monitor (mug/m^3)
#'     	\item sodium_ion. Sodium ion PM2.5 measured at monitor (mug/m^3)
#'     	\item OC. Organic carbon PM2.5 measured at monitor (mug/m^3)
#'     	\item nitrate. Nitrate PM2.5 measured at monitor (mug/m^3)
#'     	\item elemental_carbon. Elemental carbon PM2.5 measured at monitor (mug/m^3)
#'     	\item sulfate. Sulfate PM2.5 measured at monitor (mug/m^3)
#'}
#'
#' @docType data
#' @keywords datasets
#' @format A list of length 85 (named by monitor; see monitorsNE.rda)containing PM2.5 constituent concentrations for each of 85 monitors.  Each element is a data.frame with N rows and 26 variables, where the rows correspond to the dates of data for each monitor.  
#' @name consConc
NULL
