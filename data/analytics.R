library('affy')   # Affymetrix pre-processing
library('limma') 

unzip('E-MEXP-850.raw.1.zip');
candidategenes = read.table('genenames.csv', sep=',', header=T);
names=c("array01.txt","array02.txt","array03.txt","array04.txt","array05.txt",
        "array06.txt", "array07.txt","array08.txt","array09.txt","array10.txt",
        "array11.txt","array12.txt","array13.txt","array14.txt", "array15.txt", 
        "array16.txt");


metadata = read.table('E-MEXP-850.sdrf.txt', header=TRUE, sep='\t');

RG<-read.maimages(files=names, columns=list(Rf="ArrayWoRx:S_Cy5", 
                                            Gf="ArrayWoRx:S_Cy3", Rb="ArrayWoRx:B_Cy5",
                                            Gb="ArrayWoRx:B_Cy3"));
# show(RG);
# summary(RG$R)
# boxplot(data.frame(log2(RG$Gb)),main="Green background");
# imageplot(log2(RG$Gb[,1]),RG$printer);
names(RG$printer)
# 
# # Add array list information
probes <- read.table("Probes.txt",header=TRUE,sep="\t",as.is=TRUE);
spottypes <- readSpotTypes("spottypes.txt")
# 
RG$genes <- probes;
RG$printer <- getLayout(probes);
RG$genes$Status <- controlStatus(spottypes,RG);

# Normalise
MA <- normalizeWithinArrays(RG);

# fit <- lmFit(MA);
# fit <- eBayes(fit);
# topTable(fit)$Name;
# names(MA)
# MA$M[1:10]
#j = which(MA$genes$Name == 'G01_o232_plate_16')
#i= !is.na(match(MA$genes$Name, candidategenes$Probe_name));

#Search forcandidate genes
i = c();
for(name in candidategenes$Probe_name)
{
  j = which(MA$genes$Name == name);
  i = append(i,j)
  #print(paste(name,'_', j));
}

# Create dataset with candidate genes
m = cbind(MA$genes[i,],M=MA$M[i,],A=MA$A[i,]);

# Change array colnames
anames = c();
for(l in 8:39)
{
 
  anames = append(anames, strsplit(colnames(m[l]),'[.]')[[1]][2]);
}

colnames(m)[8:39] = anames;
colnames(candidategenes)
colnames(m)
m = merge(candidategenes[,1:3], m, by.x='Probe_name', by.y='Name');

timepoint = unique(metadata$Factor.Value..time.);
# Average arrays per reps
for(t in timepoint)
{
  wl = c();
  i = which(metadata$Factor.Value..time. == t);
  filename = print(metadata$FileName[i]);
  for(w in filename)
  {
    wl = append(wl,strsplit(w,'[.]')[[1]][1])
  }
  o = match(as.character(wl), colnames(m));
  m = data.frame(m, tempt= apply(m[,o], 1, mean));
  l = which(colnames(m) == 'tempt');
  colnames(m)[l] = paste('t_',t);
  
}

