## jtk: not necessary -- scripts should be run in clean R process (R CMD BATCH or similar)
## rm(list=ls())
library("affy")   # Affymetrix pre-processing
library("limma") 
#library(som);


softwarebin = function()
{
  
  break();
  # Analise single channel
  tname = c();
  db = c();
  for(fname in names)
  {
    dset = read.table(fname, header = TRUE, sep = "\t");
    r = mean(dset$ArrayWoRx.S_Cy5)/mean(dset$ArrayWoRx.S_Cy3); # ratio cy5/cy3 (red/green)
    
    dcy5 = sapply(dset$ArrayWoRx.S_Cy5, function(x) x/r); # Normalise using ratio
    i = intersect(which(metadata$FileName == fname), 
                  which(metadata$Parameter.Value..Label.used..1 == "Cy5")); # Search arrays 
    tname = append(tname, as.vector(metadata$Factor.Value..time.[i]));
    
    dcy3 = sapply(dset$ArrayWoRx.S_Cy3, function(x) x/r);
    i = intersect(which(metadata$FileName == fname), 
                  which(metadata$Parameter.Value..Label.used..1 == "Cy3"));
    
    tname = append(tname, as.vector(metadata$Factor.Value..time.[i]));
    db = cbind(db, dcy5, dcy3);
  }
  
  
  source("https://bioconductor.org/biocLite.R")
  biocLite("affy")
  #
  colnames(db) = tname; # chance colnames to timepoints
  db = MeanPerCol(timepoint, db); # Get mean per timepopint, reduce to one column
  d = data.frame(dset$Reporter.identifier, db); ## Add probe name
  colnames(d)[2:9] = as.vector(timepoint);# Change colnames to single timepoint
  # End dataframe assembly
  
  probes = data.frame(probes, IDNAME = sapply(probes$ID, 
                                              function(x) paste("ebi.ac.uk:MIAMExpress:Reporter:A-MEXP-556.",x,sep = "" )));
  
  candidategenes = candidategenes[which(candidategenes$type == 1),]; # Search for candidate genes out of paper list
  # Search for candidate genes
  i = c();
  for(name in candidategenes$Probe_name)
  {
    j = match(name, probes$Name);
    i = append(i,j)
  }
  
  z = match(probes$IDNAME[i], d$dset.Reporter.identifier); # Search for gene names 
  EmpiricalCandidateGenes = cbind(probe_name = probes$Name[i], d[z,]); # Create a dataset
  EmpiricalCandidateGenes = merge(candidategenes, EmpiricalCandidateGenes, 
                                  by.x = "Probe_name", by.y = "probe_name");
  
  write.table(EmpiricalCandidateGenes, file = "empirical_candidate.csv", sep = ",", row.names = FALSE); # Export data
  
}



# dEMO
MeanPerCol <- function(timepoint_list, data)
{
  copydata <- data;
  db1 <- c();
  for(tp in timepoint_list)
  {
    i <- which(colnames(copydata) == tp);
    ds <- apply(copydata[,i], 1, mean);
    db1 <- cbind(db1, ds);
  }
  return(db1);
}

readEmpiricalData <- function(fname)
{
  return(read.csv(fname, stringsAsFactors = FALSE, check.names = FALSE));
}

plotProfile <- function(e, probeName, geneName, xname, yname)
{
  p <- as.numeric(e)
  barplot(p, main = geneName, names.arg = names(e), las = 2, col = "red", xlab = xname, ylab = yname );
}


readData <- function(filenameList, probesFname, spottypesFname)
{
  columnList <- list(Rf = "ArrayWoRx:S_Cy5", 
                     Gf = "ArrayWoRx:S_Cy3",
                     Rb = "ArrayWoRx:B_Cy5",
                     Gb = "ArrayWoRx:B_Cy3");
  RG <- read.maimages(files = filenameList, columns = columnList);
  probes <- read.table(probesFname, header = TRUE, sep = "\t", as.is = TRUE);
  spottypes <- readSpotTypes(spottypesFname);
  RG$genes <- probes;
  RG$printer <- getLayout(probes);
  RG$genes$Status <- controlStatus(spottypes, RG);
  RG$weights<-as.matrix(read.table(file="weights.txt",header=TRUE,sep="\t"));
  
  return(RG);
}


readGregersenData <- function(genetableFname, metadataFname, probesFname, spottypesFname)
{
  g <- list();
  ## FIXME: should do this only when necessary?
  unzip("E-MEXP-850.raw.1.zip");
  g$genetable <- read.table(genetableFname, sep = ",", header = TRUE, stringsAsFactors = FALSE);
  g$metadata <- read.table(metadataFname, header = TRUE, sep = "\t", stringsAsFactors = FALSE);
  filenameList <- unique(g$metadata$FileName);
  g$RG <- readData(filenameList, probesFname, spottypesFname);
  return(g);
}

  
normaliseData <- function(data)
{
  RG <- data;
  #MA <- normalizeWithinArrays(RG);
  print('hehere1')
  RGnb <- backgroundCorrect(RG,method="none")
  print('hehere2')
  MAnorm <- normalizeWithinArrays(RGnb,span=0.4)
  print('hehere3')
  MAnormxscale<- normalizeBetweenArrays(MAnorm,method="scale")
  print('hehere4')
  MAgenesonly<-MAnormxscale[MAnormxscale$genes$Status=="Gene",]
  print('hehere5')
  #MA <- normalizeWithinArrays(RG, bc.method="none");
  
  return(MAgenesonly);
  
}


changeNameArray = function(targettable)
{
  
  for(cn in colnames(targettable))
  {
    #for(n in lc)
    #{
      for(i in 1:8)
      {
        #o = which(levels(targettable[[cn]]) == n)
        levels(targettable[[cn]])[i] = paste('t',i, sep='');
      }
        
    #}
  }
  return(targettable)
}

getTargets <- function(metadata)
{
  
  targetList <- list()
  for(dye in unique(metadata$Parameter.Value..Label.used..1))
  {
    for(an in unique(metadata$FileName))
    {
      i <- intersect(which(metadata$Parameter.Value..Label.used..1 == dye),
                    which(metadata$FileName == an));
      
      targetList[[dye]][an] <- as.character(metadata$Factor.Value..time.[i])
      
    }    
  }
  return(data.frame(targetList))
}

  
fitModel <- function(MA, design)
{
  fit <- lmFit(MA, design);
  fit <- eBayes(fit);
  return(fit);
}


# create model and get data
createLM <- function(g, tref)
{
  MA <- normaliseData(g$RG);
  targets <- getTargets(g$metadata);
  targets = changeNameArray(targets)
  print(targets)
  ## FIXME: is -5dap the appropriate reference?
  design <- modelMatrix(targets, ref = tref)
  f <- fitModel(MA, design);
  return(f);
}


# plot genes
plotGenes <- function(fit, genetable, probeList)
{
  for(probeNameProfile in probeList$probeName)
  {
    i <- which(fit$genes$Name == probeNameProfile);
    j <- which(genetable$Probe_name == probeNameProfile);
    geneDescription <- genetable$Category_description_blasthit[j];
    #png('plot.png')
    plotProfile(fit$coefficients[i,], probeNameProfile, geneDescription, "timepoint", "intensity (log10)");
#dev.off()
    if(readline(sprintf("%s %s", probeNameProfile, geneDescription)) == "q" ) { break();}
  }
}


## jtk: following http://stackoverflow.com/questions/3548090/facet-grid-problem-input-string-1-is-invalid-in-this-locale
## bioconductor seems a bit locale-sensitive (yes, locales are ugly stuff anyway)
Sys.setlocale(locale="C")

g <- readGregersenData("genenames.csv", "E-MEXP-850.sdrf.txt", "Probes.txt", "spottypes.txt");

fit <- createLM(g, 't1') # fit using t1 as reference
topTable(fit)$Name;
candidategenes <- c("G01_o232_plate_16");
probeList <- read.csv("probelist.csv");
plotGenes(fit, g$genetable, probeList);

# Gregersen
eb <-eBayes(fit);
cont.matrix2 <- makeContrasts(t8vst2=t8-t2,levels=fit$design)  #I only tested between the very late stage and the flag leaf around pollination.
fit1c <- eBayes(contrasts.fit(fit, cont.matrix2))  # this orders the genes according to their differential expression between t8 and t2 (late and early).
plotGenes(fit1c, g$genetable, probeList);
