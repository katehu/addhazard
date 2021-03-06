#' An hypothetical two-phase sampling dataset based on \link{nwtsco} dataset from
#' the National Wilms Tumor Study (NWTS)
#' @format A data frame with 3915 rows and 15 variables:
#' \describe{
#' We create a hypothetical two-phase sampling (stratified sampling) dataset.
#' In this dataset, only the subjects who have relapses and some of the controls have
#' their samples sent to the central laboratory for more accurate histology examination.
#'
#' Institutional histology is examined in the local hospital. It is usually  available for
#' all the samples. Central histology is more expensive to obtain since the samples have to be
#' sent to the central laboratory and the work requires experienced lab scientists.
#' It is reasonable to assume only some samples were tested for central
#' histology. Following the two-phase sampling design, we selected subjects for central histology
#' measurements based on their outcomes and institutional histology results.
#'
#' \item{trel}{Time to relapse or last date seen (yr), continuous}
#' \item{tsur}{Time to death or last date seen (yr), continuous}
#' \item{relaps}{Indicator of relapse,
#'                0 = Alive no prior relapse when last seen,
#'                1 = Relapsed after trel years}
#' \item{dead}{Indicator of death,
#'               0 = Alive when last seen,
#'               1 = Died after tsur years}
#' \item{study}{NWTS study,
#'              3 = NWTS-3,
#'              4 = NWTS-4}
#'\item{stage}{Stage of disease,
#'             1=I,
#'             2=II,
#'             3=III,
#'             4=IV}
#' \item{histol}{Central Path histology,
#'              0 = Favorable (FH) and the subject is selected into the phase II subsample (in.ph2 = 1),
#'              1 = Unfavorable (UH) and the subject is selected into the phase II subsample (in.ph2 = 1),
#'              NA = subject is NOT selected into the phase II subsample (in.ph2 = 1) }
#' \item{instit}{Institutional histology,
#'               0 = Favorable (FH),
#'              1 = Unfavorable (UH)}
#' \item{age}{Age at diagnosis (yr), continuous}
#' \item{yr}{Year of diagnosis, calendar year}
#' \item{specwgt}{Weight of tumor bearing specimen, in grams (continuous)}
#' \item{tumdiam}{Diameter of tumor, in centimeters (continuous)}
#' \item{strt}{Strata,
#'              1 = Unfavorable Institutional histology and no relapse,
#'              2 = favorable Institutional histology and no relapse,
#'              3 = relapse}
#' \item{Pi}{Sampling probability,
#'              0.5 for strata =1,
#'              0.9 for strata = 2,
#'              0.9 for strata = 3}
#' \item{wts}{weight assigned to each individual, inverse of the selection probability}
#' \item{in.ph2}{Phase II membership,
#'              1 = selected into the phase II subsample,
#'              0 = not selected into the phase II subsample}
#'   }
#' @source
#'  This dataset was created based on \link{nwtsco} dataset from the National Wilms Tumor Study (NWTS)
#' @docType data
#' @keywords datasets
#' @name nwts2ph
"nwts2ph"
