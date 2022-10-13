# Motif enrichment analysis using Homer

Before running the motif enrichment analysis, the DMR files need to be convered to the propper format for homer.

In this case I'll be using a .bed file for input into the findMotifsGenome.pl pipeline

.bed file column order:
1. chromosome
2. start
3. end
4. chromsome:start-end
5. strand

-size 376 = median size of DMR region. used to randomly parse genome into 376bp segments to use as background.
-cpg = normalize CpG instead of GC content
-noweight = disable GC normalization
-p 20 = multithreading


```{R}
library(IRanges)
liver.dmr<-data.table::fread("picta_liver_dmrs_all.csv")
muscle.dmr<-data.table::fread("picta_muscle_dmrs_all.csv")

sex.dmr<-findOverlapPairs(as(muscle.sex,"GRanges"),as(liver.sex,"GRanges"))
sex.dmr<-punion(sex.dmr,fill.gap=TRUE)
sex.dmr<-as.data.frame(sex.dmr)
sex.dmr$faidx<-paste(sex.dmr$seqnames,sex.dmr$start,sep=":")
sex.dmr$faidx<-paste(sex.dmr$faidx,sex.dmr$end,sep="-")
sex.dmr<-sex.dmr[!duplicated(sex.dmr$faidx),]

#Generate background regions by concatenating liver and muscle DMRs keeping only the largest regions

auto.dmr<-findOverlapPairs(as(muscle.dmr,"GRanges"),as(liver.dmr,"GRanges"))

temp<-as.data.frame(auto.dmr)
temp$muscle<-paste(temp$first.X.seqnames,temp$first.X.start,sep=".")
temp$liver<-paste(temp$second.X.seqnames,temp$second.X.start,sep=".")
muscle.dmr$muscle<-paste(muscle.dmr$chr,muscle.dmr$start,sep=".")
liver.dmr$liver<-paste(liver.dmr$chr,liver.dmr$start,sep=".")

muscle.auto<-muscle.dmr[!(muscle.dmr$muscle %in% temp$muscle),]
liver.auto<-liver.dmr[!(liver.dmr$muscle %in% temp$liver),]
muscle.auto$faidx<-paste(muscle.auto$chr,muscle.auto$start,sep=":")
muscle.auto$faidx<-paste(muscle.auto$faidx,muscle.auto$end,sep="-")

auto.dmr<-punion(auto.dmr,fill.gap=TRUE)
auto.dmr<-as.data.frame(auto.dmr)
auto.dmr$faidx<-paste(auto.dmr$seqnames,auto.dmr$start,sep=":")
auto.dmr$faidx<-paste(auto.dmr$faidx,auto.dmr$end,sep="-")
auto.dmr<-auto.dmr[!duplicated(auto.dmr$faidx),]
temp<-muscle.auto[,c(1,2,3,4,11,10)]
colnames(temp)<-c("seqnames","start","end","width","strand","faidx")
temp2<-rbind(auto.dmr,temp)
write.table(temp2[,c(1:3,6,5)],"auto_dmr_all_homer.bed",row.names=FALSE,col.names=FALSE,quote=FALSE,sep="\t")


```
```bash
perl findMotifsGenome.pl sex_dmr_overlap_homer.bed picta_genome.fasta ./ -bg auto_dmr_all_homer.bed  -p 40

```
