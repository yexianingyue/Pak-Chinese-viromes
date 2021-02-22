#=========================================
#             fig2
#=========================================
#------------------------
# fig 2a (ok)
# Virus detection（means：mapped least 20 reads）
#------------------------
library(vegan)
library(ggplot2)
get_pan_data=function(dd, nm=nm){
  dd.curve=specaccum(dd, method = "random")
  dd.curve.data=data.frame(Sites=dd.curve$sites, Richness=dd.curve$richness, SD=dd.curve$sd)
  dd.curve.data$label=rep(nm, nrow(dd.curve.data))
  dd.curve.data
}

dt = read.table("../00.data/dna.mapping.reads",sep="\t", header=T,row.names=1, check.names = F)
map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)

dt[dt<20] = 0
dt[dt>=20] = 1

dt = dt[,-61]

pa = map[which(map$group1=="Pakistan" & map$group2=="adult"),1]
pc = map[which(map$group1=="Pakistan" & map$group2=="child"),1]
ca = map[which(map$group1=="China" & map$group2=="adult"),1]
cc = map[which(map$group1=="China" & map$group2=="child"),1]

dpa = dt[,pa]
dpc = dt[,pc]
dca = dt[,ca]
dcc = dt[,cc]


ddpa = get_pan_data(t(dpa),"pa")
ddpc = get_pan_data(t(dpc),"pc")
ddca = get_pan_data(t(dca),"ca")
ddcc = get_pan_data(t(dcc),"cc")
data = rbind(ddpa, ddpc, ddca,ddcc)


ggplot(data, aes(x=Sites, y=Richness, group=label, color=label))+
  geom_line(size=1)+
  geom_errorbar(aes(ymax = Richness + SD, ymin = Richness - SD), width = 0.25)+
  scale_color_brewer(palette = "Set3")+theme_bw()+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())



#------------------------
# fig 2d
# related reads ,
#------------------------
rm(list=ls())
library(ggplot2)
library(ggpubr)
library(reshape2)
library(vegan)
dt = read.table("../00.data/dna.mapping.reads",sep="\t", header=T,row.names=1, check.names = F)
sample_map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)

sample_map$group3 = paste(sample_map$group1, sample_map$group2,sep=" ")
dt = dt[,-61]
dt = t(t(dt)/colSums(dt))
dt = sqrt(dt)
dist = as.matrix(vegdist(t(dt), method='bray'))
mds = monoMDS(dist)
point = mds$points

data = merge(point, sample_map, by.x='row.names', by.y='sample')

ggscatter(data, x = "MDS1", y = "MDS2",
          color = "group3", palette = "npg",
          ellipse = TRUE, ellipse.type="confidence",ellipse.level = 0.95,
          mean.point = F, star.plot = TRUE, 
          ggtheme = theme_bw())


#------------------------
# fig 2f
# reads related
#------------------------
rm(list=ls())
library(ggplot2)
library(reshape2)
dt = read.table("../00.data/dna.mapping.reads",sep="\t", header=T,row.names=1, check.names = F)
sample_map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)
votu_taxo =  read.table("../00.data/votu.group.txt", sep="\t", header=T, check.names = F)[,c(1,3)]

sample_map$group3 = paste(sample_map$group1, sample_map$group2,sep=" ")
dt = dt[,-61]

dt = t(t(dt)/colSums(dt))

dm = merge(votu_taxo, dt, by.x='name', by.y='row.names')[,-1]
dm = aggregate(.~family, dm, sum)

dl = melt(dm)
data = merge(dl, sample_map, by.x='variable', by.y='sample')
ggplot(data=data, aes(x=variable, y=value, fill=family))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  facet_grid(.~group3, scale='free')+
  theme(axis.text.x = element_blank())




#==================================================
#             fig3
#==================================================

#------------------------
# fig 3a 
# 把reads map到注释到kegg的gene上面，算了一下相对丰度
#------------------------
rm(list=ls())
library(ggplot2)
library(reshape2)
dt = read.table("../00.data/dna.kegg.reads.ko.levelB.desc",sep="\t", header=T,row.names=1, check.names = F)
sample_map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)
sample_map$group3 = paste(sample_map$group1, sample_map$group2,sep=" ")
dt = dt[,-61]

dt = as.data.frame(t(t(dt)/colSums(dt)))
dt = dt[order(rowSums(dt), decreasing = T),]
dt$desc = rownames(dt)
dt$desc[11:51] = "other"
ds = aggregate(.~desc, dt, sum)
dl = melt(ds)
data = merge(dl, sample_map, by.x='variable', by.y='sample')
data$desc = factor(data$desc, level=unique(dt$desc))
ggplot(data=data, aes(x=variable, y=value, fill=desc))+
  geom_bar(stat='identity')+
  theme_bw()+
  scale_fill_brewer(palette = "Set3")+
  facet_grid(.~group3, scale='free')+
  theme(axis.text.x = element_blank())


#------------------------
# fig 3b
# 计算kegg相对丰度的clr丰度，然后比较
#------------------------

rm(list=ls())
library(ggplot2)
library(compositions)
library(vegan)

dt = read.table("../00.data/dna.kegg.reads.ko.levelB.desc",sep="\t", header=T,row.names=1, check.names = F)
sample_map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)

sample_map$group3 = paste(sample_map$group1, sample_map$group2,sep=" ")
dt = dt[,-61]
sample_map = sample_map[which(sample_map$group2=='adult'),]
dt = dt[,sample_map$sample]
dt = as.data.frame(t(t(dt)/colSums(dt)))
dt = dt[order(rowSums(dt), decreasing = T),]
dt = as.data.frame(clr(dt))

dt$desc = rownames(dt)
ds = aggregate(.~desc, dt, sum)
dl = melt(ds)
data = merge(dl, sample_map, by.x='variable', by.y='sample')
data$desc = factor(data$desc, level=unique(dt$desc))




nspecies = nrow(dt)
names = rownames(dt)
grps = unique(sample_map$group1)
result = rbind()
for(n in 1:nspecies){
  grp1 = sample_map[which(sample_map$group1 == "Pakistan"),1]
  grp2 = sample_map[which(sample_map$group1 == "China"),1]
  t1 = as.numeric(dt[n,grp1])
  t2 = as.numeric(dt[n,grp2])
  m1=mean(t1)
  m2=mean(t2)
  c1=sum(t1!=0)
  c2=sum(t2!=0)
  p=wilcox.test(t1,t2)$p.value
  temp = data.frame(n=names[n],pak_mean=m1, ch_mean=m2,pvalue=p, c1_num=c1, c2_num=c2)
  result = rbind(result, temp)
}
ss = result[which(result$pvalue < 0.05),1]
#write.table(result, "fig3b.wilcox.csv",sep=",",row.names=F)
#write.table(fdrtool(read.table("clipboard",sep="\t", header=F)[,1], statistic = 'pvalue')$qval, "clipboard",sep="\t",row.names=F)
data = data[which(data$desc %in% ss),]

ggplot(data=data, aes(x=value, y=desc, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid=element_blank())+
  scale_fill_manual(values=c("#4dbbd5", "#3c5488"))


#------------------------
# fig 3c 
# 计算kegg相对丰度的clr丰度，然后比较
#------------------------

rm(list=ls())
library(fdrtool)
library(ggplot2)
library(compositions)
library(vegan)

dt = read.table("../00.data/dna.kegg.reads.ko.levelB.desc",sep="\t", header=T,row.names=1, check.names = F)
sample_map = read.table("../00.data/sample.group",sep="\t", header=T, check.names = F)

sample_map$group3 = paste(sample_map$group1, sample_map$group2,sep=" ")
dt = dt[,-61]
sample_map = sample_map[which(sample_map$group2=='child'),]
dt = dt[,sample_map$sample]
dt = as.data.frame(t(t(dt)/colSums(dt)))
#dt = dt[order(rowSums(dt), decreasing = T),]
dt = as.data.frame(clr(dt))

dt$desc = rownames(dt)
ds = aggregate(.~desc, dt, sum)
dl = melt(ds)
data = merge(dl, sample_map, by.x='variable', by.y='sample')
data$desc = factor(data$desc, level=unique(dt$desc))




nspecies = nrow(dt)
names = rownames(dt)
grps = unique(sample_map$group1)
result = rbind()
for(n in 1:nspecies){
  grp1 = sample_map[which(sample_map$group1 == "Pakistan"),1]
  grp2 = sample_map[which(sample_map$group1 == "China"),1]
  t1 = as.numeric(dt[n,grp1])
  t2 = as.numeric(dt[n,grp2])
  m1=mean(t1)
  m2=mean(t2)
  c1=sum(t1!=0)
  c2=sum(t2!=0)
  p=wilcox.test(t1,t2)$p.value
  temp = data.frame(n=names[n],pak_mean=m1, ch_mean=m2,pvalue=p, c1_num=c1, c2_num=c2)
  result = rbind(result, temp)
}
ss = result[which(result$pvalue < 0.05),1]
#write.table(result, "fig3c.wilcox.csv",sep=",",row.names=F)
#write.table(fdrtool(read.table("clipboard",sep="\t", header=F)[,1], statistic = 'pvalue')$qval, "clipboard",sep="\t",row.names=F)
data = data[which(data$desc %in% ss),]


ggplot(data=data, aes(x=value, y=desc, fill=group1))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank())+
  scale_fill_manual(values=c("#e64b35", "#00a087"))

#------------------------
# fig 3d
# 把reads map到注释到cazy的gene上面，算了一下相对丰度
#------------------------

rm(list=ls())
library(ggplot2)
library(reshape2)
library(ggsci)
library(patchwork)
dt = read.table("../00.data/dna.cazy.uniq.f.taxo.num",sep="\t", header=F,check.names = F)
dt_sum = aggregate(.~V2, dt[,c(1,2)], sum)
dt_order = dt_sum[order(dt_sum$V1, decreasing = F),1]
dt$V2 = factor(dt$V2, level=dt_order)
ggplot(data=dt, aes(x=V1, y=V2, fill=V3))+
  geom_bar(stat='identity', position='stack')+
  theme_bw()+
  scale_fill_brewer(palette='Set2')+
  theme(panel.grid=element_blank())

dt_pie = aggregate(.~V3, dt[,c(1,3)], sum)
dt_pie$percent = round(dt_pie$V1/sum(dt_pie$V1)*100, digits = 2)
dt_pie$lab = paste(dt_pie$V3, "(",dt_pie$percent,"%)", sep="")
ggplot(data=dt_pie, aes(x='1',y=V1, fill=V3))+
  geom_bar(stat='identity', position='stack')+
  theme_bw()+
  scale_fill_brewer(palette='Set2')+
  theme(panel.grid=element_blank(),
        axis.text.x=element_blank())+
  coord_polar(theta = 'y')+
  geom_text(aes(y=dt_pie$V1/2 + c(0,cumsum(dt_pie$V1)[-length(dt_pie$V1)]),
                x=1.5,
                label=rev(lab)))


#=========================================
#             figs1
#=========================================

# 密度曲线
rm(list=ls())
dt = read.table("../00.data/figs1.data",sep="\t",check.names = F, header=T)
map = read.table("../00.data/votu.group.txt", sep="\t", header=T, check.names = F)
dm = merge(dt, map, on='name')
library(ggplot2)

library(ggpubr)
library(ggsci)
library(aplot)
dm$depth = log10(dm$depth)
dm$length = log10(dm$length+150)

p0 = ggplot(data=dm)+
  geom_point(aes(x=length,y=depth, shape=circle, color=family, size=circle))+
  theme_bw()+
  scale_shape_manual(values=c(19,15))+
  scale_color_brewer(palette = "Set3")+
  scale_size_manual(values=c(3,1.5))+
  theme(panel.grid =element_blank())

#密度
xp = ggplot(dm,aes(x=length))+
  geom_density()+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.grid =element_blank())

yp = ggplot(dm,aes(x=depth))+
  geom_density()+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal()+
  theme(axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),panel.grid =element_blank())+
  coord_flip()

p0%>%
  insert_top(xp,height = 0.1)%>%
  insert_right(yp,0.1)




#==================================================
#             fig4
#==================================================
#------------------------
# fig 4f (ok)
# 已经抽平化的RNA reads， 物种名称只选择相似度>=95的
#------------------------
rm(list=ls())
dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t",check.names = F, header=T, row.names=1)
sample_map = read.table("../00.data/sample.group", sep="\t", header=T, check.names=F)
rna_taxo = read.table("../00.data/rna.group.txt", sep="\t", header=T, check.names = F)[,c(1,4)]

sample_map$group3 = paste(sample_map$group1, sample_map$group2, sep=" ")

dm = merge(rna_taxo, dt, by.x='name',by.y='row.names', all.y=T)[,-1]
dm[is.na(dm$species),1]='unknown'
dms = aggregate( . ~ species, dm,sum)
dms = dms[order(rowSums(dms[,2:61]),decreasing = T),]
rownames(dms)=dms$species
dms=dms[,-1]
data = as.data.frame(t(t(dms)/colSums(dms)))
data$name = rownames(data)
data$name[11:28] = 'other'
data = aggregate(.~name, data, sum)
dl = melt(data)
dm = merge(dl, sample_map, by.x='variable', by.y='rna_sample')

ggplot(data=dm, aes(x=variable,y=value,fill=name))+
  geom_bar(stat='identity')+
  theme_bw()+
  theme(panel.grid = element_blank(),
        axis.text.x=element_blank())+
  facet_grid(.~group3, scale="free")+
  scale_fill_brewer(palette='Set3')

#------------------------
# fig 4g (ok)
# 已经抽平化的RNA reads， 物种名称只选择相似度>=95的
#------------------------
rm(list=ls())

#处理数据------
dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t",check.names = F, header=T, row.names=1)
sample_map = read.table("../00.data/sample.group", sep="\t", header=T, check.names=F)
rna_taxo = read.table("../00.data/rna.group.txt", sep="\t", header=T, check.names = F)[,c(1,4)]

sample_map$group3 = paste(sample_map$group1, sample_map$group2, sep=" ")
sample_map = sample_map[which(sample_map$group2=='adult'),]

dm = merge(rna_taxo, dt, by.x='name',by.y='row.names', all.y=T)[,-1]
dm[is.na(dm$species),1]='unknown'
dms = aggregate( . ~ species, dm,sum)
dms = dms[order(rowSums(dms[,2:61]),decreasing = T),]
rownames(dms)=dms$species
dms=dms[,-1]

nspecies = nrow(dms)
names = rownames(dms)
grps = unique(sample_map$rna_sample)
result = rbind()
for(n in 1:nspecies){
  grp1 = sample_map[which(sample_map$group1 == "Pakistan"),4]
  grp2 = sample_map[which(sample_map$group1 == "China"),4]
  t1 = as.numeric(dms[n,grp1])
  t2 = as.numeric(dms[n,grp2])
  m1=mean(t1)
  m2=mean(t2)
  c1=sum(t1!=0)
  c2=sum(t2!=0)
  p=wilcox.test(t1,t2)$p.value
  temp = data.frame(n=names[n],pak_mean=m1, ch_mean=m2,pvalue=p, c1_num=c1, c2_num=c2)
  result = rbind(result, temp)
}

#write.table(result, "fig4g.wilcox.csv",sep=",",row.names=F)
#write.table(fdrtool(read.table("clipboard",sep="\t", header=F)[,1], statistic = 'pvalue')$qval, "clipboard",sep="\t",row.names=F)

#画图------
sel = read.table("../00.data/fig4g.list.txt",sep="\t", header=F, check.names=F)
data = dms[sel$V1,]
data$name = row.names(data)
dl = melt(data)
dm = merge(dl, sample_map, by.x='variable', by.y='rna_sample')
dm$name = factor(dm$name, sel$V1)
ggplot(data=dm, aes(x=name, y=log10(value), fill=group1))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=c("#4dbbd5", "#3c5488"))+
  facet_grid(.~name,scale='free')






#------------------------
# fig 4h (ok)
# 已经抽平化的RNA reads， 物种名称只选择相似度>=95的
#------------------------
rm(list=ls())

#处理数据------
dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t",check.names = F, header=T, row.names=1)
sample_map = read.table("../00.data/sample.group", sep="\t", header=T, check.names=F)
rna_taxo = read.table("../00.data/rna.group.txt", sep="\t", header=T, check.names = F)[,c(1,4)]

sample_map$group3 = paste(sample_map$group1, sample_map$group2, sep=" ")
sample_map = sample_map[which(sample_map$group2=='child'),]

dm = merge(rna_taxo, dt, by.x='name',by.y='row.names', all.y=T)[,-1]
dm[is.na(dm$species),1]='unknown'
dms = aggregate( . ~ species, dm,sum)
dms = dms[order(rowSums(dms[,2:61]),decreasing = T),]
rownames(dms)=dms$species
dms=dms[,-1]

nspecies = nrow(dms)
names = rownames(dms)
grps = unique(sample_map$rna_sample)
result = rbind()
for(n in 1:nspecies){
  grp1 = sample_map[which(sample_map$group1 == "Pakistan"),4]
  grp2 = sample_map[which(sample_map$group1 == "China"),4]
  t1 = as.numeric(dms[n,grp1])
  t2 = as.numeric(dms[n,grp2])
  m1=mean(t1)
  m2=mean(t2)
  c1=sum(t1!=0)
  c2=sum(t2!=0)
  p=wilcox.test(t1,t2)$p.value
  temp = data.frame(n=names[n],pak_mean=m1, ch_mean=m2,pvalue=p, c1_num=c1, c2_num=c2)
  result = rbind(result, temp)
}

#write.table(result, "fig4h.wilcox.csv",sep=",",row.names=F)
#write.table(fdrtool(read.table("clipboard",sep="\t", header=F)[,1], statistic = 'pvalue')$qval, "clipboard",sep="\t",row.names=F)

#画图------
sel = read.table("../00.data/fig4h.list.txt",sep="\t", header=F, check.names=F)
data = dms[sel$V1,]
data$name = row.names(data)
dl = melt(data)
dm = merge(dl, sample_map, by.x='variable', by.y='rna_sample')
dm$name = factor(dm$name, sel$V1)
ggplot(data=dm, aes(x=name, y=log10(value), fill=group1))+
  geom_boxplot()+
  theme_bw()+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle=45, hjust=1))+
  scale_fill_manual(values=c("#e64b35", "#00a087"))+
  facet_grid(.~name,scale='free')





#==================================================
#             fig5
#==================================================
#------------------------
#     fig 5a(ok)
#------------------------
rm(list=ls())
library(vegan)
library(ggplot2)
dna_dt = read.table("../00.data/dna.mapping.reads",sep="\t", header=T, row.names=1,check.names=F)
rna_dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t", header=T, row.names=1, check.names=F)

dna_dt = dna_dt[,-61]
dna_dt = as.data.frame(t(t(dna_dt)/colSums(dna_dt)))
colnames(rna_dt) = gsub("R","D", colnames(rna_dt))

dna_shannon = as.data.frame(diversity(t(dna_dt), index = 'shannon'))
rna_shannon = as.data.frame(diversity(t(rna_dt), index = 'shannon'))
data = merge(dna_shannon, rna_shannon,by='row.names')
colnames(data) = c("sample","dna","rna")
ggplot(data=data, aes(x=rna,y=dna))+
  geom_point()+
  theme_bw()+
  theme(panel.grid = element_blank())

#------------------------
#     fig 5b(ok)
#------------------------
rm(list=ls())
library(vegan)
library(ade4)
dna_dt = read.table("../00.data/dna.mapping.reads",sep="\t", header=T, row.names=1,check.names=F)
rna_dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t", header=T, row.names=1, check.names=F)
sample_map = read.table("../00.data/sample.group", sep="\t", header=T, check.names=F)

sample_map$group3 = paste(sample_map$group1, sample_map$group2, sep=" ")
dna_dt = dna_dt[,-61]
dna_dt = as.data.frame(t(t(dna_dt)/colSums(dna_dt)))
rna_dt = as.data.frame(t(t(rna_dt)/colSums(rna_dt)))
colnames(rna_dt) = gsub("R","D", colnames(rna_dt))# 替换

dna_dist = as.data.frame(as.matrix(vegdist(t(dna_dt), method = 'bray')))
rna_dist = as.data.frame(as.matrix(vegdist(t(rna_dt), method = 'bray')))

dna = as.data.frame(monoMDS(dna_dist)$point)
rna = as.data.frame(monoMDS(rna_dist)$point)

dna = dna[rownames(rna),]

pro_test = protest(X = dna, Y = rna, permutations = 999)
procu = procuste(dna,rna)

D = as.data.frame(procu$tabX); D$sample = rownames(D)
R = as.data.frame(procu$tabY); R$sample=rownames(R)


dna$sample = rownames(dna)
rna$sample = rownames(rna)

data = rbind(D,  R)
data$method = rep(c("DNA","RNA"), each=60)

data = merge(data, sample_map,by='sample')



ggplot(data, aes(MDS1, MDS2, color = group3, shape = method)) + 
  geom_point(size = 3) +
  geom_line(aes(group = sample, color = group3), alpha=0.5) + 
  ggtitle("Procustes correlation = 0.3706, p=0.001") + 
  scale_color_npg() + 
  theme_bw()+ 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



#------------------------
#     fig 5c(ok)
#------------------------
rm(list=ls())
rna_dt = read.table("../00.data/rna_mapped_reads_sub44944",sep="\t", header=T, row.names=1, check.names=F)
