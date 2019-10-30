library(gdata)
library(ggplot2)
library(dplyr)
library(igraph)
library(reshape2)
require(gridExtra)
#Calibrate pixels


plot_graph <- function(exp_mat){
  
  g <- graph_from_adjacency_matrix(as.matrix(exp_mat[1:nrow(exp_mat), 
                                                     2:ncol(exp_mat)]) )
  
  return(g)
}  


sum_intensities <- function(data)
{
  
  copydata = data
  copydata = data.frame(copydata, sumint = sapply(rownames(copydata), 
              function(x) sum(copydata[x,3:131])))
  return(copydata)
  
}



convert2hex <- function(data1, data2)
{
  copydata1 = data1
  copydata2 = data2
  hexlist = c()
  
  for(i in 1:length(copydata1))
  {
    idx = which(copydata2$i == copydata1[i])
    hexlist = append(hexlist, as.character(copydata2[idx,5]))
    
  }
  return(hexlist)
}



get_subdata = function(sdata, conver_table, intens)
{

  
  #' I us want to use the same intensities # Vey crude way but it 
  #' seems they're the most relevants.
  
  u = colnames(sdata) %in% intens # Search for intens
  scc1= data.frame(idtag = sdata$idtag, 
                   DAS = sdata$DAS, sdata[,u]) # Create frame with sel intensi 
  freq = colnames(scc1)[3:dim(scc1)[2]]
  freq = as.numeric(sub("X", "", freq))
  colnames(scc1)[3:dim(scc1)[2]] = freq # Delete X in columns
  #write.table(scc1, file='t.csv', sep=',', row.names=F)
  
  mm = as.matrix(scc1[, 3:dim(scc1)[2]]) # Deletes -4
  vector_data = unmatrix(mm,byrow=T) # Convert matrix into vector
  vector_data = as.numeric(vector_data) # Convert values to numeric
  idtag = as.character(unique(scc1$idtag)) # Unique idtags
  DAS = as.character(unique(scc1$DAS)) # unique time points
  hexlist = convert2hex(freq, conver_table) #
  hexlist[c(1:11)] <- c("#C08040","#C08041","#C08042","#C08043","#C08044",
                        "#C08045","#C08046","#C08047","#C08048","#C08049",
                        "#C08040")
  hexlist[c(12:18, 20)] <- c("#004001", "#004002","#004003", "#004004",
                             "#004005","#004006", "#004007", "#004008")
  conveg = expand.grid(idtag=idtag, intensity=hexlist,DAS=DAS)
  cdata = data.frame(conveg, value=vector_data)
  
  
  
  #cdatacal = calibrate_data(cdata)
  return(cdata)
}

# Take row means

takeMeans <- function(x) 
{
  x1 = rapply(x, mean)
  x1['time'] = names(x) 
  return(x1)
}

# Plot the thing
print_plot <-  function(sdata, br, lb, o){
 
  
  p <- ggplot(sdata,aes(x = DAS, y = value, fill=intensity, group=intensity)) +
    geom_area() + geom_line(aes(y=value), position="stack") +
    scale_fill_manual(values = as.character(unique(sdata$intensity))) + 
    ylab(label = "Colour frequency") + xlab(label = "Developmental stage") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.text.y = element_text(size = 12),
          axis.title.x = element_text(size = 12),
          axis.title.y = element_text(size = 12),
          panel.background = element_rect(fill = "white", colour = "grey50"),
          legend.position="none") +
    labs(title = 'B', face ='bold') +
    #scale_x_continuous(breaks=seq(min(sdata$DAS), max(sdata$DAS), by=5), 
    #                   labels = seq(min(sdata$DAS), max(sdata$DAS), by=5)) +
    scale_x_continuous(breaks = br, labels = lb) +
    #ggtitle(unique(sdata$idtag)) + 
    geom_vline(xintercept = c(as.numeric(as.character(o[[1]])), 
                              as.numeric(as.character(o[[2]])),
                              as.numeric(as.character(o[[3]]))),
               #               as.numeric(as.character(o[[4]]))), 
               linetype="dotted", 
               color = c("red", "green", "blue"), size=1.5)
    
  print(p)
}

# subset of data
average_similar_measures = function(data, featurename)
{
  copydata = data
  copydata[[featurename]] = factor(copydata[[featurename]])
  #print(plantid)
  # Arrange new dataset
  l <- split(copydata, copydata[[featurename]])
  h <- do.call(rbind, lapply(1:length(l), function(x) takeMeans(l[x])))
  h = h[,c(-1,-2)]
  h = h[,c(dim(h)[2], 1:(dim(h)[2]-1))]
  h <- data.frame(idtag=unique(copydata$idtag), h)
  colnames(h)[3:(dim(h)[2])] <- colnames(scsub)[3:(dim(scsub)[2])]
  return(h)
}

calibratedata = function(sdata, DAS, date_list)
{
  
  i <- which(sdata[['DAS']] <= DAS[match(as.Date("2015-04-01", '%Y-%m-%d'), 
                                     date_list)])
  
  sdata[i, 3:131] = sdata[i, 3:131]*0.27
  sdata[-i, 3:131] = sdata[-i, 3:131]*0.56
  return(sdata)
}


convert2numeric = function(data, col_list)
{
  
  copydata=data
  for(cname in col_list)
  {
    copydata[[cname]] = as.numeric(as.character(copydata[[cname]]))
  }
  return(copydata)
}


get_DAS <- function(sdata, DAS, date_list){
  
  
  DAS_doodle <- rep(0, nrow(sdata))
  for(ts in unique(sdata[['timestamp']])){
    
    DAS_doodle[as.numeric(rownames(subset(sdata, timestamp == ts)))] = 
      DAS[match(as.Date(ts, format = '%Y-%m-%d'), 
                date_list)]
    
  }
  
  sdata <- cbind(sdata, DAS = DAS_doodle)
  
  return(sdata)
}


#' Plot colour image
plot_color_image <- function(idtagname_list, doodle, groundtruth, it){
  
  
  for(idtagname in idtagname_list){
    
    scsub = subset(doodle, idtag == idtagname)
    h <- subset(aggregate(scsub, by=list(scsub[['DAS']]), 
                          FUN=mean), select = -c(idtag, DAS, timestamp))
    h <- cbind(idtag = idtagname, h)
    colnames(h)[2] <- 'DAS'
    cdatacal <- get_subdata(h, it, intens)
    cdatacal$DAS <- as.numeric(as.character(cdatacal$DAS))
    o <- groundtruth[match(idtagname, groundtruth$barcode), 
                     c('gs39_DAS', 'gs55_DAS', 
                       'GS65_DAS', 'FlagLeafSenescence_DAS')]
    if(length(which(complete.cases(as.matrix(o)[1,]) == 'FALSE')) == 0){
      if(max(as.numeric(as.matrix(o)[1,])) < max(cdatacal$DAS) & 
         min(as.numeric(as.matrix(o)[1,])) > min(cdatacal$DAS)){
        colnames(o) <- c('GS39', 'GS55', 'GS65', 'FLS')
        print(idtagname)
        
        br <-  seq(min(cdatacal$DAS), max(cdatacal$DAS), by=1)
        lb <- rep('', length(br))
        lb[as.numeric(lapply(1:4, function(x) which(br == o[[x]])))] = colnames(o)
        
        tiff(paste(idtagname, '.tiff', sep = ''), height = 1000, width = 1000,
             res = 200)
        print_plot(cdatacal, br, lb, o)
        dev.off()
      }
    }
    
    #if(readline(idtagname) == 'q') { break }
  }
}

#' Get average color data
average_color_data <- function(groundtruth, doodle){
  
  #'Get average colour data
  idtagname_list <- unique(groundtruth$barcode)
  all_data <- c()
  for(idtagname in idtagname_list){
    
    scsub = subset(doodle, idtag == idtagname)
    h <- subset(aggregate(scsub, by=list(scsub[['DAS']]), 
                          FUN=mean), select = -c(idtag, DAS, timestamp))
    h <- cbind(idtag = idtagname, h)
    
    all_data <- rbind(all_data, h)
    #if(readline(idtagname) == 'q') { break }
  }
  return(all_data)
}


plot_gene_expression_data <- function(emp_data, gene_group, color_code){
  
  DT.m1 <- melt(emp_data, id.vars = 'name',
                measure.vars = names(emp_data)[-1], variable.name = 'stage',
                value.name = 'expression_level')
  
  
  DT.m1 <- data.frame(DT.m1, group = gene_group)
  
 
  p <- ggplot(DT.m1, aes(x = stage, y = expression_level,  group = name, 
                         colour= name)) +
    geom_point(size =3, aes(shape = name)) + geom_line(size=1.2) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          axis.text.x = element_text(angle = 90, hjust = 1, size = 12),
          axis.title.y = element_text(size = 12), 
          legend.text = element_text(size=10),
          legend.position = 'none') +
    #legend.position=c(0.08, 0.78)) +
    xlab(label = "Developmental stage") + ylab(label = 'Expression level (log2)') +
    labs(title = 'A', face ='bold') + ylim(0,8) +
    scale_color_manual(values=color_code)
                        
  return(p)
 
}

plot_average_colour <- function(avgdata, groundtruth){
  
  
  colnames(avgdata)[which(colnames(avgdata) == 'Group.1') ] <- 'DAS'
  avgdata <- cbind(idtag = 'avgdata', avgdata)
  cdatacal <- get_subdata(avgdata, it, intens)
  cdatacal$DAS <- as.numeric(as.character(cdatacal$DAS))
  o <- list()
  #, 'FlagLeafSenescence_DAS')]
  o[['gs39_DAS']] <- round(mean(groundtruth$gs39_DAS, na.rm = TRUE))
  o[['gs55_DAS']] <- round(mean(groundtruth$gs55_DAS, na.rm = TRUE))
  o[['GS65_DAS']] <- round(mean(groundtruth$GS65_DAS, na.rm = TRUE))
  #o[['FlagLeafSenescence_DAS']] <- round(mean(groundtruth$FlagLeafSenescence_DAS, 
  #na.rm = TRUE))
  names(o) <- c('GS39', 'GS55', 'GS65')#, 'FLS')
  
  br <-  seq(min(cdatacal$DAS), max(cdatacal$DAS), by=1)
  lb <- rep('', length(br))
  lb[as.numeric(lapply(1:length(o), function(x) which(br == o[[x]])))] = names(o)
 
  
  print_plot(cdatacal, br, lb, o)
  
}


plot_gene_graph <- function(cl, path){
  
  label_name <- c('A', 'B', 'C', 'D', 'E', 'F', 'G')
  
 
  for(i in 4:10){
    
    senes_net <- read.csv(paste(path, 'senes_t_', i, '.csv', sep=''), 
                          header = TRUE)
    o <- plot_graph(senes_net)
    l <- layout_with_sugiyama(o, layers = NULL)
   
    plot(o, vertex.size= 5*igraph::degree(o), edge.color = 'black', 
         edge.width = 1, edge.arrow.size=0.3, vertex.color = cl,  
         vertex.label.color="black",vertex.label.cex=0.9,
         vertex.label = V(o)$name,
         main = colnames(emp_data)[i+1], layout=l$layout)
    legend("topright", label_name[i-3], bty='n')
   
  }
  #dev.off()  
}
  
#' Process raw data
#' @param sdata datase
#' @param b aggregate by
process_rawdata <- function(sdata, b){
  
  
  f <- as.formula(paste('. ~', b, sep = ''))
  sdata <- aggregate(f, data = sdata, FUN= "mean")
  
  return(sdata)
  
}

#' Process raw data
#' @param sdata dataset
#' 
process_rawdata_age <- function(sdata){
  
  p <- data.frame()
  for(log2 in unique(sdata[['log2']])){
      
    p <- rbind(p, data.frame(log2 = log2, gene  = 'gene', name = 'NAM', 
                             t(sapply(4:ncol(sdata), function(x)
      mean(sdata[intersect(grep(log2, sdata[['log2']]), 
                          grep('NAM', sdata[['name']])), x])))))
  }
  
  names(p) <- names(sdata)
  sdata <- rbind(sdata, p)
  sdata <- sdata[-grep('NAM_', sdata[['name']]),]
  sdata[['name']] <- drop.levels(sdata[['name']])
  sdata[['name']] <- factor(sdata[['name']], 
                            levels = unique(sdata[['name']]))
  

  sdata <- melt(sdata, id.vars = c("log2","gene","name"),
                measure.vars = names(sdata)[-1:-3], variable.name = 'stage',
                value.name = 'expression_level')
  return(sdata)
  
}

#' Plost anthesis data
#' @param sdata dataset
#'
plot_anthesis <- function(sdata){
  
  p <- ggplot(sdata,aes(x = stage, y = expression_level, color = name)) +
    geom_point(size =3, aes(shape = name)) + geom_line(size=1.2) +
    theme(panel.background = element_rect(fill = "white", colour = "grey50"),
          axis.text.x = element_text(angle = 0, hjust = 1, size = 12),
          axis.title.y = element_text(size = 12), 
          legend.text = element_text(size=10)) +
    xlab(label = "Days after anthesis") + ylab(label = 'Expression level (log2)') +
    #labs(title = 'A', face ='bold') + ylim(0,8) +
    scale_color_manual(values=color_code) +
    labs(title = 'B', face ='bold') + ylim(0,8)
  
  return(p)

}


doodle <- read.table('doodle.csv', header=T, sep=',')
freqs <- colnames(doodle)[3:dim(doodle)[2]]
freqs <- as.numeric(sub("X", "", freqs))
colnames(doodle)[3:dim(doodle)[2]] = freqs
groundtruth <- read.table('groundtruth.csv', header=TRUE, sep=',', 
                         na.strings = c('?', 'NA'))

start_day <- '2014-10-20'
end_day <- '2015-04-30'
date_list <- d <- seq.Date(as.Date(start_day, format = '%Y-%m-%d'), 
                           as.Date(end_day, format = '%Y-%m-%d'), by='day') 

DAS <- sapply(d, function(x) x-d[1])
doodle <- get_DAS(doodle,  DAS, date_list)
doodle <- calibratedata(doodle, DAS, date_list)

# Read table with frequencies
it <- read.table('contable.csv', header=T, sep=',')
it <- data.frame(it, hex=rgb(it[,2:4], max=255))
# Select intensities
intens = c(1,4,5,16,20,21,24,25,26,27,30,31,32,50,54,55,60,108,112,116)
# Select subdataset

idx = which(groundtruth$unknown != '')
idtagname_list = unique(groundtruth$barcode[idx]) 
# 

# Plot colour plots
#plot_color_image(idtagname_list, doodle, groundtruth, it)
all_data <- average_color_data(groundtruth, doodle)
write.table(all_data, file='all_calib_data.csv', sep=',', quote = F, 
            row.names=F)

gene_group <- c('Chlorophyll','vrn-B3','WRKY','MYB', 'NAM')
# Plot gene expression data
emp_data <- read.csv('data/empiricaldata_v1.csv', header=TRUE, skip = 1)
emp_data <- process_rawdata(emp_data[, 3:13], 'name') 



color_code <- c("#999999", "#E69F00", "#56B4E9", "#009E73", 
      "#D55E00") #, "#b20019","#F0E442", 



break
tiff('senes_net.tiff', width=900, height = 800, res = 150)
par(mfrow = c(3,3), mar=c(0, 1.5, 1, 0))
plot_gene_graph(color_code, 'c:/Users/x992029/share/')
dev.off()
# Plot
avgdata <- subset(aggregate(all_data, by=list(all_data[['Group.1']]), 
                            FUN=mean), select = -c(idtag, Group.1))
tiff(paste('average', '.tiff', sep = ''), height = 1000, width = 1000,
     res = 200)
plot_average_colour(avgdata, groundtruth)
dev.off()
## Plot gene 


p1 <- plot_gene_expression_data(emp_data, gene_group, color_code)

anthesis_data <- read.csv('data/empiricaldata_postanthesis.csv', header = TRUE)
anthesis_data <- process_rawdata_age(anthesis_data)
anthesis_data$stage <- as.numeric(gsub('X', '', anthesis_data$stage))

p2 <- plot_anthesis(subset(anthesis_data, log2 == 'tpm'))


tiff('expression_prof.tiff', height = 800, width = 1500, res = 180)
grid.arrange(p1, p2, ncol=2)
dev.off()

