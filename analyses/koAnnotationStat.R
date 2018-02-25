library(ggplot2)
library(gridExtra)
setwd("/Users/trvinh/work/OLD/koAnnotation/toolValidation")

##### DISTRIBUTION PLOT FOR T_FAS_KO #####
koRefIn <- read.csv("FAS/koRefspec/koThreshold_mean.list.OLD", quote='', sep='\t', header=F)
koRefMean = round(mean(koRefIn$V2),3)

p <- ggplot(koRefIn, aes(x=V2)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=koRefIn,mapping=aes(xintercept=mean(V2)),linetype="dashed", size=1, color="red")+
  geom_text(data=koRefIn,aes(x=koRefMean-0.25, y=2500, label=paste("mean = ",koRefMean)), size=8, vjust=-0.4, color="red")
p <- p + theme_minimal() + 
          theme(legend.position = "top",
               legend.title = element_blank(),
               legend.text = element_text(size=15),
               axis.text = element_text(size=20),
               axis.title = element_text(size=20)) +
          labs(x="KO specific FAS cutoff",y="count")
p

### 2 examples for low and high T_FAS KO groups
scaleFUN0 <- function(x) sprintf("%.0f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)

ko_K00542 <- read.csv("FAS/koRefspec/K00542.FAS", quote='', sep='\t', header=F)
mean_K00542 <- round(mean(ko_K00542$V4),3)
p1 <- ggplot(ko_K00542, aes(x=V4)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=ko_K00542,mapping=aes(xintercept=mean_K00542),linetype="dashed", size=1, color="red")+
  geom_text(data=ko_K00542,aes(x=mean_K00542, y=8, label=paste("mean = ",mean_K00542)), size=6,angle=90, vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p1 <- p1 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of K00542",y="count")
p1

ko_K07888 <- read.csv("FAS/koRefspec/K07888.FAS", quote='', sep='\t', header=F)
str(ko_K07888)
mean_K07888 <- round(mean(ko_K07888$V4),3)
p2 <- ggplot(ko_K07888, aes(x=V4)) +
  geom_histogram(binwidth=.003, alpha=.5, position="identity") +
  geom_vline(data=ko_K07888,mapping=aes(xintercept=mean_K07888),linetype="dashed", size=1, color="red")+
  geom_text(data=ko_K07888,aes(x=mean_K07888, y=6, label=paste("mean = ",mean_K07888)), size=6, angle=90,vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0)
p2 <- p2 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of K07888",y="count")
p2

grid.arrange(p1,p2,ncol=2)


##### ANNOTATED YEAST PROTEINS STAT #####

##### FAS score distribution of supported and unsupported orthologs ##### 
scaleFUN0 <- function(x) sprintf("%.0f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)

allAnno <- read.csv("sacce_rbh_checkco/anno/sacce_anno.list.NEW.KO.fasList", quote='', sep='\t',header=F)
mean_allAnno <- round(mean(allAnno$V1),3)
p <- ggplot(allAnno, aes(x=V1)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=allAnno,mapping=aes(xintercept=mean_allAnno),linetype="dashed", size=1, color="red")+
  geom_text(data=allAnno,aes(x=mean_allAnno-0.1, y=1500, label=paste("mean = ",mean_allAnno)), size=4, vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p <- p + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of HamFAS orthologs",y="count")
# p

supported <- read.csv("sacce_rbh_checkco/anno/sacce_anno.list.NEW.KO.withINPA.FAS", quote='', sep='\t',header=F)
mean_supported <- round(mean(supported$V1),3)
p1 <- ggplot(supported, aes(x=V1)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=supported,mapping=aes(xintercept=mean_supported),linetype="dashed", size=1, color="red")+
  geom_text(data=supported,aes(x=mean_supported-0.1, y=1500, label=paste("mean = ",mean_supported)), size=4,vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p1 <- p1 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of supported orthologs",y="count")
# p1

unsupported <- read.csv("sacce_rbh_checkco/anno/sacce_anno.list.NEW.KO.noINPA.FAS", quote='', sep='\t',header=F)
mean_unsupported <- round(mean(unsupported$V1),3)
p2 <- ggplot(unsupported, aes(x=V1)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=unsupported,mapping=aes(xintercept=mean_unsupported),linetype="dashed", size=1, color="red")+
  geom_text(data=unsupported,aes(x=mean_unsupported-0.1, y=100, label=paste("mean = ",mean_unsupported)), size=4,vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p2 <- p2 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of unsupported orthologs",y="count")
# p2

grid.arrange(p,p1,p2,ncol=1)

##### ANNOTATED YEAST PROTEINS STAT #####

##### FAS score distribution of supported and unsupported orthologs ##### 
scaleFUN0 <- function(x) sprintf("%.0f", x)
scaleFUN1 <- function(x) sprintf("%.1f", x)

allUnk <- read.csv("sacce_rbh_checkco/un_anno/sacce_unknown.list.NEW.KO.fasList", quote='', sep='\t',header=F)
mean_allUnk <- round(mean(allUnk$V1),3)
p1 <- ggplot(allUnk, aes(x=V1)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=allUnk,mapping=aes(xintercept=mean_allUnk),linetype="dashed", size=1, color="red")+
  geom_text(data=allUnk,aes(x=mean_allUnk-0.1, y=80, label=paste("mean = ",mean_allUnk)), size=4, vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p1 <- p1 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS scores of HamFAS orthologs",y="count")
# p1

hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/sacce_unknown.list.NEW.KO.fasList_hamfasOnly", quote='', sep='\t',header=F)
mean_hamfasOnly <- round(mean(hamfasOnly$V1),3)
p2 <- ggplot(hamfasOnly, aes(x=V1)) +
  geom_histogram(binwidth=.01, alpha=.5, position="identity") +
  geom_vline(data=hamfasOnly,mapping=aes(xintercept=mean_hamfasOnly),linetype="dashed", size=1, color="red")+
  geom_text(data=hamfasOnly,aes(x=mean_hamfasOnly-0.1, y=40, label=paste("mean = ",mean_hamfasOnly)), size=4, vjust=-0.4, color="red") +
  scale_y_continuous(labels=scaleFUN0) +
  scale_x_continuous(labels=scaleFUN1)
p2 <- p2 + theme_minimal() + 
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.text = element_text(size=15),
        axis.text = element_text(size=15),
        axis.title = element_text(size=20)) +
  labs(x="FAS score of HamFAS-only orthologs",y="count")
# p2

grid.arrange(p1,p2,ncol=1)

##### compare length and # of domains with non-HamFAS-only #####
hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/hamfasOnly/sacce_unknown.list.NEW.KO.hamfasOnly.LengthDomain", quote='', sep='\t',header=T)
non_hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/hamfasOnly/sacce_unknown.list.NEW.KO.others.LengthDomain", quote='', sep='\t',header=T)
hamfasOnly$type <- "HamFas_only"
non_hamfasOnly$type <- "Others"

df <- rbind(hamfasOnly,non_hamfasOnly)
df[df$length==0,]
df[df$length == min(df$length),]

p_len <- ggplot(df, aes(length, colour=type))
p_len <- p_len + geom_density(alpha=0.55)
p_len <- p_len + theme_minimal() +
        labs(x="sequence length")
p_len

p_domain <- ggplot(df, aes(domains, colour=type))
p_domain <- p_domain + geom_density(alpha=0.55)
p_domain <- p_domain + theme_minimal() +
  labs(x="# of PFAM domains")
p_domain

##### analyze HamFAS-only #####
hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/hamfasOnly/sacce_unknown.list.NEW.KO.hamfasOnly.LengthDomain", quote='', sep='\t',header=T)
non_hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/hamfasOnly/sacce_unknown.list.NEW.KO.others.LengthDomain", quote='', sep='\t',header=T)
hamfasOnly$type <- "HamFas_only"
non_hamfasOnly$type <- "Others"

head(hamfasOnly)
df <- rbind(hamfasOnly,non_hamfasOnly)
df[df$length==0,]
df[df$length == min(df$length),]

dfSub <- df[df$origin == "bacteria" | df$origin == "archaea",]
dfSub

dfSub[seqID == "sacce:4062",]

p_len <- ggplot(dfSub, aes(length, colour=type))
p_len <- p_len + geom_density(alpha=0.55)
p_len <- p_len + theme_minimal() +
  labs(x="sequence length")
p_len

p_domain <- ggplot(dfSub, aes(domains, colour=type))
p_domain <- p_domain + geom_density(alpha=0.55)
p_domain <- p_domain + theme_minimal() +
  labs(x="# of PFAM domains")
p_domain

##### analyze origin of annotations #####
library(plyr) # for count function
unkAnno <- read.csv("sacce_rbh_checkco/un_anno/sacce_unknown.list.NEW.KO.LengthDomain", quote='', sep='\t',header=T)
anno <- read.csv("sacce_rbh_checkco/anno/sacce_anno.list.NEW.KO.LengthDomain", quote='', sep='\t',header=T)
hamfasOnly <- read.csv("sacce_rbh_checkco/un_anno/hamfasOnly/sacce_unknown.list.NEW.KO.hamfasOnly.LengthDomain", quote='', sep='\t',header=T)

count(unkAnno$origin)
count(anno$origin)
count(hamfasOnly$origin)

##### number of orthologs #####
orthoAnno <- read.csv("sacce_rbh_checkco/anno/sacce_anno.list.NEW.KO.orthoCount", quote='', sep='\t',header=F)
orthoUnk <- read.csv("sacce_rbh_checkco/un_anno/sacce_unknown.list.NEW.KO.orthoCount", quote='', sep='\t',header=F)
orthoHamfas <- read.csv("sacce_rbh_checkco/un_anno/sacce_unknown.list.NEW.KO.hamfasOnly.orthoCount", quote='', sep='\t',header=F)

orthoAnno$type <- "annotated"
orthoUnk$type <- "unannotated"
orthoHamfas$type <- "unannotated_hamfasOnly"

nrow(orthoHamfas)
nrow(orthoAnno[orthoAnno$V2>1,])
nrow(orthoUnk[orthoUnk$V2>1,])
nrow(orthoHamfas[orthoHamfas$V2>1,])

orthoCountDf <- rbind(orthoAnno,orthoUnk)
orthoCountDf <- rbind(orthoCountDf,orthoHamfas)

p_count <- ggplot(orthoCountDf, aes(V2, colour=type))
p_count <- p_count + geom_density(alpha=0.55)
p_count <- p_count + theme_minimal() +
  labs(x="Number of orthologs")
p_count


##### analyze PPI (connectivity) #####
connectAnno <- read.csv("sacce_rbh_checkco/anno/connectivity/sacce_anno.list.NEW.KO.PPI", quote='', sep='\t',header=T)
connectUnk <- read.csv("sacce_rbh_checkco/un_anno/connectivity/sacce_unknown.list.NEW.KO.PPI", quote='', sep='\t',header=T)
connectUnkHamfas <- read.csv("sacce_rbh_checkco/un_anno/connectivity/sacce_unknown.list.NEW.KO.hamfasOnly.PPI", quote='', sep='\t',header=T)

connectAnno$type <- "annotated"
mean(connectAnno$nodeDegree)
connectUnk$type <- "unannotated"
mean(connectUnk$nodeDegree)
connectUnkHamfas$type <- "unannotated_hamfasOnly"
mean(connectUnkHamfas$nodeDegree)

nrow(connectUnk[connectUnk$nodeDegree > 10,])
nrow(connectUnk)

connectDf <- rbind(connectAnno,connectUnk)
connectDf <- rbind(connectDf,connectUnkHamfas)
head(connectDf)
connectDf[connectDf$nodeDegree <= 100,]

p_connect <- ggplot(connectDf, aes(nodeDegree, colour=type))
p_connect <- p_connect + geom_density(alpha=0.55)
p_connect <- p_connect + theme_minimal() +
  labs(x="PPI degree")
p_connect

p_connect2 <- ggplot(connectDf, aes(nodeDegree,..density.., colour = type)) + geom_freqpoly(binwidth = 50)
p_connect2 <- p_connect2 + theme_minimal() +
  labs(x="PPI degree")
p_connect2

# c1 <- ggplot(connectAnno, aes(nodeDegree)) +geom_histogram(binwidth=.5)
# c2 <- ggplot(connectUnk, aes(nodeDegree)) +geom_histogram(binwidth=.5)
# c3 <- ggplot(connectUnkHamfas, aes(nodeDegree)) +geom_histogram(binwidth=.5)
# grid.arrange(c1,c2,c3,ncol=3)

##### analyze pathway (connectivity) #####
pathAnno <- read.csv("sacce_rbh_checkco/anno/connectivity/sacce_anno.list.NEW.KO.ko2path", quote='', sep='\t',header=T)
pathUnk <- read.csv("sacce_rbh_checkco/un_anno/connectivity/sacce_unknown.list.NEW.KO.ko2path", quote='', sep='\t',header=T)
pathUnkHamfas <- read.csv("sacce_rbh_checkco/un_anno/connectivity/sacce_unknown.list.NEW.KO.hamfasOnly.ko2path", quote='', sep='\t',header=T)

pathAnno$type <- "annotated"
pathUnk$type <- "unannotated"
pathUnkHamfas$type <- "unannotated_hamfasOnly"

nrow(pathAnno[pathAnno$nb_pathways == 0,])
nrow(pathAnno[pathAnno$nb_pathways >= 1,])
nrow(pathAnno)

nrow(pathUnk[pathUnk$nb_pathways == 0,])
nrow(pathUnk[pathUnk$nb_pathways >= 1,])
nrow(pathUnk)

nrow(pathUnkHamfas[pathUnkHamfas$nb_pathways == 0,])
nrow(pathUnkHamfas[pathUnkHamfas$nb_pathways >= 1,])
nrow(pathUnkHamfas)


pathDf <- rbind(pathAnno,pathUnk)
pathDf <- rbind(pathDf,pathUnkHamfas)
pathDf <- pathDf[,c("KO","nb_pathways","type")]
pathDf

p_path <- ggplot(pathDf, aes(nb_pathways, colour=type))
p_path <- p_path + geom_density(alpha=0.5)
p_path <- p_path + theme_minimal() +
  labs(x="number of patwhays")
p_path

p_path2 <- ggplot(pathDf, aes(nb_pathways,..density.., colour = type)) + geom_freqpoly(binwidth = 1)
p_path2 <- p_path2 + theme_minimal() +
  labs(x="number of patwhays")
p_path2

# p1 <- ggplot(pathAnno, aes(nb_pathways)) +geom_histogram(binwidth=.5, colour="red")
# p1 <- ggplot(pathUnk, aes(nb_pathways)) +geom_histogram(binwidth=.5)
# p1 <- ggplot(pathUnkHamfas, aes(nb_pathways, fill=type)) +geom_histogram(binwidth=.5)
# grid.arrange(p1,p2,p3,ncol=3)
