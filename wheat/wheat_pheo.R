makeDayNames <- function(dayList)
{
  return(sprintf("day%03d", dayList));
}


transexprTable <- function(fname, numTimesteps)
{
  transexprCmd <- sprintf("transexpr -n %d %s", numTimesteps, fname);
  d <- read.table(pipe(transexprCmd), header = TRUE, sep = "");
  d <- d[, c(1, grep("\\.avg", colnames(d)))];
  colnames(d) <- sub("\\.avg", "", colnames(d));
  return(d);
}


findTransexprProfile <-function(tt, factorName, dayList = NULL)
{
  if (is.null(dayList))
  {
    dayList <- 1L:nrow(tt);
  }
  p <- as.numeric(tt[[factorName]][dayList]);
  names(p) <- makeDayNames(dayList);
  return(p);
}


findEmpiricalProfile <- function(e, probeName, dayList)
{
  b <- rownames(e) == probeName;
  if (sum(b) != 1)
  {
    stop(sprintf("probe %s nonexistent or ambiguous", probeName));
  }
  dayNames <- makeDayNames(dayList);
  p <- as.numeric(e[b, dayNames]);
  names(p) <- dayNames;
  return(p);
}


readEmpiricalTable <- function(fname, dayList, probeNameColumn = "probes.Name.i.")
{
  e <- read.csv(fname, stringsAsFactors = FALSE);
  ## FIXME: should check that dayList and number of columns match up
  rownames(e) <- e[[probeNameColumn]];
  colnames(e)[3:10] <- makeDayNames(dayList);
  return(e);
}


readWheat <- function(dayList)
{
  numDays <- max(dayList);
  w <- list();
  w$e <- readEmpiricalTable("data/csv/empiric_candidate_v1.csv", dayList);
  w$t <- transexprTable("wheat_pheo.trl", numDays);
  return(w)
}


plotProfile <- function(p, ...)
{
  barplot(p, las = 2, ...);
}


plotProfiles <- function(w, probeFactorTable, dayList, waitFunc = function(m) { readline(m); })
{
  opar <- par(no.readonly = TRUE);
  par(mfrow = c(2, 1));
  probeNameList <- rownames(probeFactorTable);
  for (probeName in probeNameList)
  {
    factorName <- probeFactorTable[probeName, "factorName"];
    pEmpirical <- findEmpiricalProfile(w$e, probeName, dayList);
    pTransexpr <- findTransexprProfile(w$t, factorName, dayList);
    plotProfile(pEmpirical);
    plotProfile(pTransexpr);
    waitFunc(sprintf("%s:%s", probeName, factorName));
  }
  par(opar);
}


# tentative mapping from days past (before) pollination to days of life based on
# http://www.extension.umn.edu/agriculture/small-grains/growth-and-development/spring-wheat/ , fig. 1
dayList <- c(58L, 63L, 73L, 83L, 88L, 93L, 99L, 105L);
probeFactorTable <- data.frame(
  factorName = c("pheophorbide_oxygenase"),
  stringsAsFactors = FALSE
  );
rownames(probeFactorTable) <-  c("TA-PaO_clone_39");

w <- readWheat(dayList);
plotProfiles(w, probeFactorTable, dayList);
