rm(list=ls())
library('affy')   # Affymetrix pre-processing
library('limma') 
#library(som);


softwarebin = function()
{
  
  break();
  # Analise single channel
  tname = c();
  db = c();
  for(fname in names)
  {
    dset = read.table(fname, header=T, sep='\t');
    r = mean(dset$ArrayWoRx.S_Cy5)/mean(dset$ArrayWoRx.S_Cy3); # ratio cy5/cy3 (red/green)
    
    dcy5 = sapply(dset$ArrayWoRx.S_Cy5, function(x) x/r); # Normalise using ratio
    i = intersect(which(metadata$FileName == fname), 
                  which(metadata$Parameter.Value..Label.used..1 == 'Cy5')); # Search arrays 
    tname = append(tname, as.vector(metadata$Factor.Value..time.[i]));
    
    dcy3 = sapply(dset$ArrayWoRx.S_Cy3, function(x) x/r);
    i = intersect(which(metadata$FileName == fname), 
                  which(metadata$Parameter.Value..Label.used..1 == 'Cy3'));
    
    tname = append(tname, as.vector(metadata$Factor.Value..time.[i]));
    db = cbind(db, dcy5, dcy3);
  }
  
  
  
  #
  colnames(db) = tname; # chance colnames to timepoints
  db = MeanPerCol(timepoint, db); # Get mean per timepopint, reduce to one column
  d = data.frame(dset$Reporter.identifier, db); ## Add probe name
  colnames(d)[2:9] = as.vector(timepoint);# Change colnames to single timepoint
  # End dataframe assembly
  
  probes = data.frame(probes, IDNAME = sapply(probes$ID, 
                                              function(x) paste('ebi.ac.uk:MIAMExpress:Reporter:A-MEXP-556.',x,sep='' )));
  
  candidategenes = candidategenes[which(candidategenes$type == 1),]; # Search for candidate genes out of paper list
  # Search for candidate genes
  i = c();
  for(name in candidategenes$Probe_name)
  {
    j = match(name, probes$Name);
    i = append(i,j)
  }
  
  z = match(probes$IDNAME[i], d$dset.Reporter.identifier); # Search for gene names 
  EmpiricalCandidateGenes = cbind(probe_name=probes$Name[i], d[z,]); # Create a dataset
  EmpiricalCandidateGenes = merge(candidategenes, EmpiricalCandidateGenes, 
                                  by.x = 'Probe_name', by.y='probe_name');
  
  write.table(EmpiricalCandidateGenes, file = 'empirical_candidate.csv', sep=',', row.names = F); # Export data
  
}



# dEMO
MeanPerCol = function(timepoint_list, data)
{
  copydata = data;
  db1 = c();
  for(tp in timepoint_list)
  {
    i = which(colnames(copydata) == tp);
    ds = apply(copydata[,i], 1, mean);
    db1 = cbind(db1, ds);
  }
  return(db1);
}

readEmpiricalData <- function(fname)
{
  return(read.csv(fname, stringsAsFactors = FALSE, check.names = FALSE));
}

plotProfile <- function(e, probeName, geneName, xname, yname)
{
  
  p = as.numeric(e)
  barplot(p, main = geneName, names.arg = names(e), las = 2, col='red', xlab = xname, ylab = yname );
}



readData = function(names, probes, spottypes)
{
  
  RG = read.maimages(files=names, columns=list(Rf="ArrayWoRx:S_Cy5", 
                                               Gf="ArrayWoRx:S_Cy3", Rb="ArrayWoRx:B_Cy5",
                                               Gb="ArrayWoRx:B_Cy3"));
  RG$genes <- probes;
  RG$printer <- getLayout(probes);
  RG$genes$Status <- controlStatus(spottypes,RG);
  return(RG);
}

  
normaliseData = function(data)
{
  RG = data;
  MA <- normalizeWithinArrays(RG);
  return(MA);
  
}


getDesignMatrix = function(timepoint)
{
  
  design = list()
  for(dye in unique(metadata$Parameter.Value..Label.used..1))
  {
    for(an in unique(metadata$FileName))
    {
      i = intersect(which(metadata$Parameter.Value..Label.used..1 == dye),
                    which(metadata$FileName == an));
      
      design[[dye]][an] = as.character(metadata$Factor.Value..time.[i])
      
    }    
  }
    
  return(data.frame(design))
}

  
fitModel = function(MA, design)
{
  
  fit <- lmFit(MA, design);
  fit <- eBayes(fit);
  return(fit);

}

# create model and get data
createLM = function(names, probes, spottypes, metadata)
{
  RG = readData(names, probes, spottypes);
  MA = normaliseData(RG);
  m = getDesignMatrix(metadata);
  design = modelMatrix(m,ref="-5")
  f = fitModel(MA, design);
  return(f);
  
}



unzip('E-MEXP-850.raw.1.zip');
genetable = read.table('genenames.csv', sep=',', header=T, stringsAsFactors = TRUE);
metadata = read.table('E-MEXP-850.sdrf.txt', header=TRUE, sep='\t', stringsAsFactors = TRUE);
probes <- read.table("Probes.txt",header=TRUE,sep="\t",as.is=TRUE);
spottypes <- readSpotTypes("spottypes.txt");


fit = createLM(unique(metadata$FileName), probes, spottypes, metadata);
topTable(fit)$Name;
candidategenes= c('G01_o232_plate_16');


# plot genes
for(probeNameProfile in candidategenes)
{
  i = which(fit$genes$Name == probeNameProfile);
  j = which(genetable$Probe_name == probeNameProfile);
  plotProfile(fit$t[i,], probeNameProfile, genetable$Category_description_blasthit[j], 'timepoint', 'intensity (log10)');
  if(readline(probeNameProfile) == 'q' ) { break();}
  
}
