library('affy')   # Affymetrix pre-processing
library('limma') 
#library(som);


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



ProcessDataLimma = function(names, probes, spottypes)
{
  
  RG = read.maimages(files=names, columns=list(Rf="ArrayWoRx:S_Cy5", 
                                               Gf="ArrayWoRx:S_Cy3", Rb="ArrayWoRx:B_Cy5",
                                               Gb="ArrayWoRx:B_Cy3"));
  show(RG);
  summary(RG$R)
  boxplot(data.frame(log2(RG$Gb)),main="Green background");
  
  # # Add array list information
  
  RG$genes <- probes;
  RG$printer <- getLayout(probes);
  RG$genes$Status <- controlStatus(spottypes,RG);
  
  # Normalise
  MA <- normalizeWithinArrays(RG);
  MA <- normalizeBetweenArrays(MA, method="Aquantile")
  
  
  # fit <- lmFit(MA);
  # fit <- eBayes(fit);
  # topTable(fit)$Name;
  # names(MA)
  # MA$M[1:10]
  #j = which(MA$genes$Name == 'G01_o232_plate_16')
  #i= !is.na(match(MA$genes$Name, candidategenes$Probe_name));
}

unzip('E-MEXP-850.raw.1.zip');
candidategenes = read.table('genenames.csv', sep=',', header=T);
names=c("array01.txt","array02.txt","array03.txt","array04.txt","array05.txt",
        "array06.txt", "array07.txt","array08.txt","array09.txt","array10.txt",
        "array11.txt","array12.txt","array13.txt","array14.txt", "array15.txt", 
        "array16.txt");

metadata = read.table('E-MEXP-850.sdrf.txt', header=TRUE, sep='\t');
timepoint = unique(metadata$Factor.Value..time.);
probes <- read.table("Probes.txt",header=TRUE,sep="\t",as.is=TRUE);
spottypes <- readSpotTypes("spottypes.txt")



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


for(probeNameProfile in EmpiricalCandidateGenes$Probe_name)
{
  i = which(EmpiricalCandidateGenes$Probe_name == probeNameProfile);
  sdata = log10(EmpiricalCandidateGenes[i, 11:18]);
  genename = as.character(EmpiricalCandidateGenes$Category_description_blasthit[i]);
  plotProfile(sdata, probeNameProfile, genename, 'timepoint', 'intensity (log10)');
  if(readline(probeNameProfile) == 'q' ) { break();}
  
}
