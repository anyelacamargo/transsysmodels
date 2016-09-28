library("ggplot2");
## FIXME: what do we need gdata for?
## can't we get rid of this? and of ggplot2 as well, preferably?
library("gdata");


readEmpirical <- function(fname, cutoff)
{
  sc <- read.table(fname, header = TRUE, sep=",");
  vc = c();
  for(i in 3:dim(sc)[2])
  {
    vc = append(vc,mean(sc[,i]))
  }
  v = which(vc <= 10);
  scc1 = sc[,c(-3, -v)];
  return(scc1);
}


cumulativePlot <- function(d)
{
  l <- list();
  for (i in 3:ncol(d))
  {
    print(i);
    ## print(rowSums(d[, 3:i]));
    if (i == 3)
    {
      l[[colnames(d)[i]]] <- d[[i]];
    }
    else
    {
      l[[colnames(d)[i]]] <-  rowSums(d[, 3:i]);
    }
  }
  colorList <- rainbow(length(l));
  x <- 1:nrow(d);
  yMax <- max(l[[length(l)]]);
  for (i in 1:length(l))
  {
    print(l[[i]]);
    if (i == 1)
    {
      print(as.character(d[["timestamp"]]));
      plot(x, l[[i]], type = "l", col = colorList[i], axes = FALSE, ylim = c(0, max(l[[length(l)]])));
      axis(1, at = 1:nrow(d), labels = as.character(d[["timestamp"]]), las = 2);
      axis(2);
      box();
    }
    else
    {
      lines(x, l[[i]], col = colorList[i]);
    }
  }
}


readPixelFreqSeries <- function(fnameList, colormapFname)
{
  a <- NULL;
  for (fname in fnameList)
  {
    cmd <- sprintf("pnmremap -mapfile=\'%s\' %s | pnmnoraw | ./pixelcount -m %s", colormapFname, fname, colormapFname);
    message(cmd);
    d <- read.table(pipe(cmd), header = TRUE, sep = "\t");
    d1 <- t(d);
    if (is.null(a))
    {
      a <- d1;
    }
    else
    {
      a <- rbind(a, d1);
    }
  }
  print(nrow(a));
  print(fnameList);
  rownames(a) <- fnameList;
  return(a);
}


### deprecated
rgb2hex <- function(maxIntensity, rgb)
{
  ## FIXME: should check for 'rgb' prefix or similar...?
  r <- as.double(substr(rgb, 4, 4)) / maxIntensity;
  g <- as.double(substr(rgb, 5, 5)) / maxIntensity;
  b <- as.double(substr(rgb, 6, 6)) / maxIntensity;
  r <- as.integer(round(r * 256));
  g <- as.integer(round(g * 256));
  b <- as.integer(round(b * 256));
  r[r > 255L] <- 255L;
  g[g > 255L] <- 255L;
  b[b > 255L] <- 255L;
  return(sprintf("#%02x%02x%02x", r, g, b));
}


getProfile <- function(cdata, idtag, intensity)
{
  b <- (cdata$idtag == idtag) & (cdata$intensity == intensity);
  timeLabel <- as.character(cdata$time[b]);
  n <- cdata$value[b];
  o <- order(timeLabel);
  n <- n[o];
  names(n) <- timeLabel[o];
  return(n);
}


allProfileList <- function(cdata, idtag)
{
  l <- list();
  for (intensity in sort(unique(cdata$intensity)))
  {
    l[[intensity]] <- getProfile(cdata, idtag, intensity);
  }
  return(l);
}


allProfileBarplot <- function(cdata, idtag, ignoreColors = character(), ...)
{
  p <- allProfileList(cdata, idtag);
  p <- p[!(names(p) %in% ignoreColors)];
  d <- t(as.data.frame(p));
  
  barplot(d, col = names(p));
}


allProfileCorrelationDistance <- function(a1, a2)
{
  d <- 0.0;
  for (color in unique(c(names(a1), names(a2))))
  {
    if (!(color %in% names(a1)) || !(color %in% names(a2)))
    {
      dCol <- 1.0;
    }
    else
    {
      p1 <- a1[[color]];
      p2 <- a2[[color]];
      if ((var(p1) == 0.0) || (var(p2) == 0.0))
      {
        if ((var(p1) == 0.0) && (var(p2) == 0.0))
        {
          dCol <- 0.0;
        }
        else
        {
          dCol <- 1.0;
        }
      }
      else
      {
        dCol <- 1.0 - cor(p1, p2);
      }
    }
    ## message(sprintf("%s: %f", color, dCol));
    d <- d + dCol;
  }
  return(d);
}


### deprecated
pfs2cdata <- function(pfs, maxIntensity)
{
  cdata <- NULL;
  intensityList <- colnames(pfs);
  for (dstep in rownames(pfs))
  {
    message(dstep);
    d <- data.frame(idtag = rep("dummyId", length(intensityList)), intensity = intensityList, time = rep(dstep, length(intensityList)), value = as.numeric(pfs[dstep, ]));
    if (is.null(cdata))
    {
      cdata <- d;
    }
    else
    {
      cdata <- rbind(cdata, d);
    }
  }
  cdata$intensity <- rgb2hex(maxIntensity, as.character(cdata$intensity));
  return(cdata);
}


## map continuous values in [0.0, 1.0] to discrete values in [0, 255]
## this implementation seems "consistent" in the sense that applying it to
## matpab_map_128.csv results in the colours found in seg.png
## Notice that a few colours in contable.csv are off by one
continuousTo8bit <- function(x)
{
  x8 <- as.integer(x * 256.0);
  x8[x8 > 255L] <- 255L;
  return(x8);
}


readCsvColormap <- function(csvFilename)
{
  return(read.csv(csvFilename, header = FALSE, stringsAsFactors = FALSE));
}


readRgbCsvPalette <- function(rgbCsvFname)
{
  d <- readCsvColormap(rgbCsvFname);
  if (ncol(d) != 3)
  {
    stop("csv needs to be in exactly 3 columns for r, g and b");
  }
  r <- continuousTo8bit(d[[1]]);
  g <- continuousTo8bit(d[[2]]);
  b <- continuousTo8bit(d[[3]]);
  return(data.frame(i = seq(along = r), hexcol = sprintf("#%02x%02x%02x", r, g, b)));
}


readContablePalette <- function(contableName)
{
  contable <- read.csv(contableName, stringsAsFactors = FALSE);
  return(data.frame(i = contable$i, hexcol = sprintf("#%02x%02x%02x", contable$R, contable$G, contable$B), stringsAsFactors = FALSE));
}


pixcolToHexcol <- function(pixcol)
{
  if (!all(substring(pixcol, 1, 3) == "pix"))
  {
    stop("pixcol contains strings not starting with \"pix\"");
  }
  return(sub("pix", "#", pixcol));
}


hexcolToPixcol <- function(hexcol)
{
  if (!all(grepl("#[0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f][0-9A-Fa-f]", hexcol)))
  {
    stop("hexcol contains invalid strings");
  }
  return(sprintf("pix%s", substring(hexcol, 2, 7)));
}


isPixColumnHeader <- function(columnHeader)
{
  return(substring(columnHeader, 1L, 3L) == "pix");
}


findPixcolList <- function(d)
{
  h <- colnames(d);
  return(h[isPixColumnHeader(h)]);
}


doodleTimeToJulian <- function(s, startDay = 0)
{
  if (is.null(s))
  {
    stop("null date string");
  }
  daysPerMonth <- c(31L, 28L, 31L, 30L, 31L, 30L, 31L, 31L, 30L, 31L, 30L, 31L);
  if (!all(substring(s, 2, 5) == "2015"))
  {
    stop("function currently works for 2015 only");
  }
  m <- as.integer(substring(s, 7, 8));
  d <- as.integer(substring(s, 10, 11));
  j <- d;
  for (i in seq(along = daysPerMonth))
  {
    b <- m > i;
    j[b] <- j[b] + daysPerMonth[i];
  }
  return(j - startDay);
}


readDoodle <- function(doodleFname, rgbCsvFname, startDay)
{
  d <- read.csv(doodleFname, header = TRUE, stringsAsFactors = FALSE);
  contable <- readRgbCsvPalette(rgbCsvFname);
  nColoursContable <- nrow(contable);
  nColoursDoodle <- ncol(d) - 2L;
  if (nColoursDoodle > nColoursContable)
  {
    warning(sprintf("truncating %d columns", nColoursDoodle - nColoursContable));
    d <- d[, 1:(ncol(d) + nColoursContable - nColoursDoodle)];
  }
  colnames(d)[3L:ncol(d)] <- hexcolToPixcol(contable$hexcol);
  j <- 
  d1 <- d[, 1:2];
  d1[["timepoint"]] <- doodleTimeToJulian(d$timestamp, startDay);
  d <- cbind(d1, d[, 3:ncol(d)]);
  return(d);
}


singlePlantDoodle <- function(d, idtag)
{
  return(d[d$idtag == idtag, ]);
}


isTimepointUnique <- function(d)
{
  return(length(unique(d$timepoint)) == nrow(d));
}


doodleCorrelationDivergence <- function(d1, d2)
{
  if (nrow(d1) != nrow(d2))
  {
    stop("incompatible doodle sets (length differ)");
  }
  if (!isTimepointUnique(d1))
  {
    stop("timepoint values in d1 not unique");
  }
  if (!isTimepointUnique(d2))
  {
    stop("timepoint values in d2 not unique");
  }
  pixList1 <- findPixcolList(d1);
  pixList2 <- findPixcolList(d2);
  if ((length(pixList1) != length(pixList2)) || (length(intersect(pixList1, pixList2)) != length(pixList1)))
  {
    stop("incompatible pixel colour lists");
  }
  if (!all(d1$timepoint == d2$timepoint))
  {
    stop("incompatible timepoint sets");
  }
  d1pix <- d1[, pixList1];
  d2pix <- d2[, pixList1];
  dcd <- 0.0;
  for (pix in pixList1)
  {
    x <- as.numeric(d1pix[[pix]]);
    y <- as.numeric(d2pix[[pix]]);
    if ((var(x) == 0.0) || (var(y) == 0.0))
    {
      if ((var(x) == 0.0) && (var(y) == 0.0))
      {
        dcd <- dcd + 0.0;
      }
      else
      {
        dcd <- dcd + 1.0;
      }
    }
    else
    {
      dcd <- dcd + 1.0 - cor(x, y);
    }
  }
  return(dcd);
}


doodleTimepointProjection <- function(d1, d2)
{
  t1 <- d1$timepoint;
  t2 <- d2$timepoint;
  if (any(!(t2 %in% t1)))
  {
    print(t1);
    print(t2);
    print(!(t2 %in% t1));
    print(t2[!(t2 %in% t1)]);
    stop("timepoints in d2 are not a subset of timepoints in d1");
  }
  b <- t1 %in% t2;
  return(d1[b, ]);
}


plotPixcolPalette <- function(pixcol)
{
  cols <- pixcolToHexcol(pixcol);
  nRows <- floor(sqrt(length(cols))) + 1;
  nCols <- ceiling(length(cols) / nRows);
  message(sprintf("nRows = %f, nCols = %f", nRows, nCols));
  xy <- expand.grid(1:nCols, 1:nRows);
  i <- seq(along = pixcol);
  x <- xy[[1]][i];
  y <- xy[[2]][i];
  plot(x, y, col = cols);
  text(x, y, pixcol, pos = 4);
}


makeTimepointCoverageTable <- function(doodle)
{
  timepointList <- sort(unique(doodle$timepoint));
  n <- rep(0L, length(timepointList));
  for (idtag in unique(doodle$idtag))
  {
    singlePlantTimepointList <- unique(singlePlantDoodle(doodle, idtag)$timepoint);
    n <- n + as.integer(timepointList %in% singlePlantTimepointList);
  }
  return(data.frame(timepoint = timepointList, numPlants = n));
}


makeNumTimepointTable <- function(doodle)
{
  idtagList <- unique(doodle$idtag);
  return(data.frame(numTimepoints = sapply(idtagList, function(idtag) { nrow(singlePlantDoodle(doodle, idtag)); })));
}


makeSinglePlantProjectionList <- function(doodle, refIdtag)
{
  refDoodle <- singlePlantDoodle(doodle, refIdtag);
  sppList <- list();
  for (idtag in unique(doodle$idtag))
  {
    s <- singlePlantDoodle(doodle, idtag);
    if (all(refDoodle$timepoint %in% s$timepoint))
    {
      sppList[[idtag]] <- doodleTimepointProjection(s, refDoodle);
    }
  }
  return(sppList);
}


correlationDivergenceMatrix <- function(doodleList)
{
  numPlants <- length(doodleList);
  m <- matrix(as.numeric(NA), nrow = numPlants, ncol = numPlants);
  rownames(m) <- names(doodleList);
  colnames(m) <- names(doodleList);
  for (idtag1 in names(doodleList))
  {
    for (idtag2 in names(doodleList))
    {
      m[idtag1, idtag2] <- doodleCorrelationDivergence(doodleList[[idtag1]], doodleList[[idtag2]]);
    }
  }
  return(m);
}


cdmDemo <- function()
{
  doodle <- readDoodle("doodle.csv", "matpab_map_128.csv", -9L);
  sppList <- makeSinglePlantProjectionList(d, "W8-001111");
  plot(hclust(as.dist(cdm)));
  return(invisible(list(doodle = doodle, sppList = sppList)));
}

### deprecated
doodle2cdata <- function(doodle, idtagName)
{
  warning("obsolete function -- use doodleToCdata");
  freqs <- colnames(doodle)[3:dim(doodle)[2]];
  freqs <- as.numeric(sub("X", "", freqs));
  colnames(doodle)[3:dim(doodle)[2]] <- freqs;
  ## idtagname <- "W8-114112";
  ix <- which(doodle$idtag == idtagName);
  scsub <- doodle[ix,];
  ## 
  it <- read.table("contable.csv", header = TRUE, sep = ",");
  it <- data.frame(it, hex = rgb(it[,2:4], max = 255));
  ## I us want to use the same intensities # Vey crude way but it seems they're the most relevants.
  intens <- c(1,4,5,16,20,21,24,25,26,27,30,31,32,50,54,55,60,108,112,116);

  u <- colnames(scsub) %in% intens;
  scc1 =  data.frame(idtag = scsub$idtag, time = scsub$timestamp, scsub[,u]);
  freq <- colnames(scc1)[3:dim(scc1)[2]];
  freq <- as.numeric(sub("X", "", freq));
  colnames(scc1)[3:dim(scc1)[2]] <- freq;
  write.table(scc1, file = "t.csv", sep = ",", row.names = FALSE);

  mm <- as.matrix(scc1[-4, 3:dim(scc1)[2]]);
  vector_data <- unmatrix(mm,byrow = TRUE);
  vector_data <- as.numeric(vector_data);
  idtag <- as.character(unique(scc1$idtag))
  time <- as.character(unique(scc1$time))

  hexlist <- c();
  for(i in 1:length(freq))
  {
    idx <- which(it$i == freq[i])
    hexlist <- append(hexlist, as.character(it[idx, 5]));
  }
  g <- expand.grid(idtag = idtag, intensity = hexlist,time = time);
  cdata <- data.frame(g, value = vector_data);
  return(cdata);
}


makeHexPalette <- function(rgbLevelList)
{
  rgbTable <- expand.grid(b = rgbLevelList, g = rgbLevelList, r = rgbLevelList);
  return(data.frame(i = 0:nrow(rgbTable), hexcol = sprintf("#%02x%02x%02x", rgbTable$r, rgbTable$g, rgbTable$b)));
}


### deprecated
doodleToCdata <- function(doodle, palette, idtagName = NULL)
{
  colorIndexDoodle <- as.integer(substring(colnames(doodle)[3:ncol(doodle)], 2L));
  colorIndex <- integer();
  colorList <- character();
  for (i in seq(along = colorIndexDoodle))
  {
    b <- palette$i == colorIndexDoodle[i];
    if (sum(b) == 0L)
    {
      warning(sprintf("no rgb color for palette index %d", colorIndexDoodle[i]));
    }
    else if (sum(b) > 1L)
    {
      stop(sprintf("found %d palette rows with i = %d", sum(b), colorIndexDoodle[i]));
    }
    else
    {
      colorList <- c(colorList, as.character(palette$rgb[i]));
      colorIndex <- c(colorIndex, i);
    }
  }
  numColors <- length(colorIndex);
  if (is.null(idtagName))
  {
    doodleSub <- doodle;
  }
  else
  {
    doodleSub <- doodle[doodle$idtag == idtagName, ];
  }
  idtagList <- character();
  intensityList <- character();
  timeList <- character();
  valueList <- integer();
  for (i in 1L:nrow(doodleSub))
  {
    valueList <- c(valueList, as.numeric(doodleSub[i, colorIndex]));
    idtagList <- c(idtagList, rep(as.character(doodleSub$idtag[i]), numColors));
    timeList <- c(timeList, rep(as.character(doodleSub$time[i]), numColors));
    intensityList <- c(intensityList, colorList);
    ## message(sprintf("lengths: valueList: %d, idtagList: %d, timeList: %d,intensityList: %d", length(valueList), length(idtagList), length(timeList), length(intensityList)));
  }
  return(data.frame(idtag = idtagList, intensity = intensityList, time = timeList, value = valueList, stringsAsFactors = FALSE));
}


doodleCdataList <- function(doodle)
{
  l <- list();
  for (idtagName in as.character(doodle$idtag))
  {
    message(idtagName);
    l[[idtagName]] <- doodle2cdata(doodle, idtagName);
  }
  return(l);
}



readLsysDoodle <- function(lsysBasename, colormapFname)
{
  fnamePatternPrefix <- sprintf("%s_d", lsysBasename);
  fnamePattern <- sprintf("%s[0-9]+\\.ppm", fnamePatternPrefix);
  # print(fnamePattern);
  fnameList <- dir(pattern = fnamePattern);
  ## FIXME: hard-coded to assume 3 digit timestamps
  timepointList <- as.integer(substring(fnameList, nchar(fnamePatternPrefix) + 1L, nchar(fnamePatternPrefix) + 3L));
  lsysPfs <- readPixelFreqSeries(fnameList, colormapFname);
  lsysDoodle <- cbind(data.frame(filename = rownames(lsysPfs), idtag = lsysBasename, timepoint = timepointList), lsysPfs);
  return(lsysDoodle);
}


### deprecated
readLsysCdata <- function(lsysBasename, colormapFname)
{
  # FIXME: hardcoded fnamePattern
  fnamePatternPrefix <- sprintf("%s_d", lsysBasename);
  fnamePattern <- sprintf("%s.+\\.ppm", fnamePatternPrefix);
  # print(fnamePattern);
  fnameList <- dir(pattern = fnamePattern);
  lsysPfs <- readPixelFreqSeries(fnameList, colormapFname);
  lsysCdata <- pfs2cdata(lsysPfs, 4);
  ## FIXME: position of time step determined by crude positional hack
  timestepStrPos <- nchar(lsysBasename) + 3;
  timeStepStr <- substring(as.character(lsysCdata$time), nchar(fnamePatternPrefix) + 1, nchar(fnameList) - 4);
  lsysCdata$time <- as.factor(sprintf("t%s", timeStepStr));
  return(lsysCdata);
}


showDistanceTable <- function(allData)
{
  for (wheatModelName in names(allData$lsysProfilesList))
  {
    message(sprintf("distance(empirical, %s) = %6.2f", wheatModelName, allProfileCorrelationDistance(allData$doodleProfiles, allData$lsysProfilesList[[wheatModelName]])));
  }
}


doodleLsys <- function()
{
  wheatModelNameList <- c("singleshoot", "alwaysgreen", "senescence");
  lsysCdataList <- list();
  lsysProfilesList <- list();
  for (wheatModelName in wheatModelNameList)
  {
    lsysCdata <- readLsysCdata(sprintf("wheat_%s", wheatModelName), "rgb4colormap.pnm");
    lsysCdataList[[wheatModelName]] <- lsysCdata;
    lsysProfiles <- allProfileList(lsysCdata, "dummyId");
    lsysProfilesList[[wheatModelName]] <- lsysProfiles;
  }
  doodleCdata <- doodle2cdata("doodle.txt");
  ## arbitrary selection of 7 time points, to match the 7 time points in the lsys data
  doodleTime <- c(" 2015-01-17", " 2015-01-23", " 2015-01-28", " 2015-02-03", " 2015-02-08", " 2015-03-05", " 2015-03-24");
  doodleTimeMap <- sort(levels(lsysCdata$time));
  names(doodleTimeMap) <- doodleTime;
  doodleCdata <- doodleCdata[doodleCdata$time %in% names(doodleTimeMap), ];
  doodleCdata$time <- doodleTimeMap[as.character(doodleCdata$time)];
  doodleCdata$intensity <- as.character(doodleCdata$intensity);
  doodleProfiles <- allProfileList(doodleCdata, "W8-114112");
  allData <- list(lsysCdataList = lsysCdataList, doodleCdata = doodleCdata, lsysProfilesList = lsysProfilesList, doodleProfiles = doodleProfiles);
  showDistanceTable(allData);
  return(invisible(allData));
}


## FIXME: unnecessary, does nothing -- delete?
findExtremes <- function(d)
{
  print(min(d))
  
}


# FIXME: canonicalise (use "<-" assignment, boolean constants, double quotes etc.)
# FIXME: separate data from code

createAverageperDate <- function(sdata)
{
  g <- aggregate(cbind(X0,X1,X2,X3,X4,X5,X6,X7,
                      X8,X9,X10,X11,X12,X13,X14,X15,X16,X17,
                      X18,X19,X20,X21,X22,X23,X24,X25,X26,X27,
                      X28,X29,X30,X31,X32,X33,X34,X35,X36,X37,
                      X38,X39,X40,X41,X42,X43,X44,X45,X46,X47,
                      X48,X49,X50,X51,X52,X53,X54,X55,X56,X57,
                      X58,X59,X60,X61,X62,X63,X64,X65,X66,X67,
                      X68,X69,X70,X71,X72,X73,X74,X75,X76,X77,
                      X78,X79,X80,X81,X82,X83,X84,X85,X86,X87,
                      X88,X89,X90,X91,X92,X93,X94,X95,X96,X97,
                      X98,X99,X100,X101,X102,X103,X104,X105,X106,X107,
                      X108,X109,X110,X111,X112,X113,X114,X115,X116,X117,
                      X118,X119,X120,X121,X122,X123,X124,X125,X126,X127) ~ timestamp, data = sdata, FUN = 'mean');
  return(g);
}


calibratedata <- function(data)
{
  copydata <- data;
  for(tdate in unique(copydata$time))
  {
    t = which(copydata$time == tdate);
    for(i in 3:(ncol(copydata) - 2))
    {
      if(as.Date(tdate) <= as.Date("2015-04-01")) 
      {
        copydata[t, i] <- copydata[t, i] * 0.27;
      }
      else
      { 
        copydata[t, i] <- copydata[t, i] * 0.46;
      }
    }
  }    
  return(copydata);
}


convert2numeric <- function(data, col_list)
{
  copydata <- data;
  for(cname in col_list)
  {
    copydata[[cname]] = as.numeric(as.character(copydata[[cname]]))
  }
  return(copydata);
}


transformData <- function(filename)
{
  cdata <- read.table(filename, header = TRUE, sep = ",");
  m <- split(cdata, cdata$timestamp);
  i <- 1;
  subdata <- m[[i]];
  s <- createAverageperDate(subdata);
  for(subdata in m[2:length(m)])
  {
    s1 <- createAverageperDate(subdata);
    s <- merge(s, s1, all = TRUE);
  }
  return(s);
}


get_subdata <- function(data, conver_table, intens)
{
  
  copydata <- data;
  # I us want to use the same intensities # Vey crude way but it seems they're the most relevants.
  
  u <- colnames(copydata) %in% intens; # Search for intens
  scc1 <- data.frame(idtag = copydata$idtag, 
                   time = copydata$time, copydata[ ,u]); # Create frame with sel intensi
  ## FIXME: dim(scc1)[2] looks like equivalent to ncol to me
  freq <- colnames(scc1)[3:ncol(scc1)];
  freq <- as.numeric(sub("X", "", freq));
  colnames(scc1)[3:ncol(scc1)] <- freq; # Delete X in columns
  mm <- as.matrix(scc1[, 3:ncol(scc1)]); # Deletes -4
  vector_data <- unmatrix(mm, byrow=T); # Convert matrix into vector
  vector_data <- as.numeric(vector_data); # Convert values to numeric
  idtag <- as.character(unique(scc1$idtag)); # Unique idtags
  time <- as.character(unique(scc1$time)) # unique time points
  hexlist <- convert2hex(freq, conver_table); # 
  conveg <- expand.grid(idtag=idtag, intensity=hexlist,time=time);
  cdata <- data.frame(conveg, value=vector_data);
  return(cdata);
}


## FIXME: what are data1, data2?
convert2hex <- function(freq_list, colormap)
{
  copydata1 <- freq_list;
  copydata2 <- colormap;
  hexlist <- c();
  
  for(i in copydata1)
  {
    hexlist <- append(hexlist, as.character(copydata2[i, 'hex']));
  }
  return(hexlist)
}


print_plot <- function(data, tname)
{
  copydata <- data;
  p <- ggplot(copydata,aes(x = time, y = value, fill = intensity, group = intensity)) +
    geom_area() + geom_line(aes(ymax=value), position="stack") +
    scale_fill_manual(values = as.character(unique(copydata$intensity))) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + ggtitle(tname);
  print(p);
}


### deprecated (can use standard rgb() function directly)
convert2col <- function(data)
{
  
  v = as.numeric(data);
  m = max(v);
  if(max(v) == 0)
  {
    m = 1;
  }
  #return(rgb(v[1], v[2], v[3],  maxColorValue=m));
  return(rgb(v[1], v[2], v[3]));
}


### deprecated
# Convert hexa to RGB
convert2RGB <- function(data)
{
  t <- c();
  for(frehex in data)
  {
    t <- rbind(t, col2rgb(frehex)[1:3]);
  }
  colnames(t) <- c('R', 'G', 'B')
  return(t);
}


### deprecated
plotDoodle <- function()
{
  system('make');
  colormap <- read.table("rgb5colormap.ppm", header = FALSE, sep = " ", skip=2);
  colormap <- data.frame(colormap, 
                         hex = sapply(rownames(colormap), 
                                      function(x) convert2col(colormap_matlab[x,1:3])));
  
  colormap_Matlab = read.table('matpab_map_128.csv', header=F, sep=',');
  
  
  colormap_Matlab <- data.frame(colormap_Matlab, 
                         hex = sapply(rownames(colormap_Matlab), 
                                      function(x) convert2col(colormap_Matlab[x,1:3])));
  
  # Convert hexa to rgb and export file
  write.table(merge(colormap_Matlab, convert2RGB(colormap_Matlab$hex), by.x='row.names', by.y='row.names'), 
              file='map.csv', sep=',');
  
                         
  ## Select intensities to avoid getting black and blue from background
  intens = c(1, 4, 5, 16, 20, 21, 24, 25, 26, 27, 30, 31, 32, 50, 54, 55, 60, 108, 112, 116);
  filename <- "doodle.csv";
  avgData <- transformData(filename);
  avgData <- data.frame('w8', avgData);
  colnames(avgData) <- c("idtag", "time", as.character(1:(ncol(avgData)-2)));
  avgData <- convert2numeric(avgData, colnames(avgData)[3:(ncol(avgData)-2)]);
  avgData <- calibratedata(avgData);
  cdatacal <- get_subdata(avgData, colormap_Matlab, intens);
  print_plot(cdatacal, "average");
  write.table(avgData, file = "avgData.csv", sep = ",", quote = FALSE, row.names = FALSE)
}


findColors <- function(colormap, colorIndexList)
{
  return(colormap[colorIndexList, ]);
}


## example usage: csvToPpm("matpab_map_128.csv", "matpab_map_128.ppm");
csvToPpm <- function(csvFilename, ppmFilename)
{
  d <- readCsvColormap(csvFilename);
  if (ncol(d) != 3)
  {
    stop("csv needs to be in exactly 3 columns, r, g and b");
  }
  r <- continuousTo8bit(d[[1]]);
  g <- continuousTo8bit(d[[2]]);
  b <- continuousTo8bit(d[[3]]);
  write(c("P3", sprintf("%d 1 255", nrow(d)), sprintf("%d %d %d", r, g, b)), file = ppmFilename);
}


#plotDoodle()
