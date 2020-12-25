
library(circlize)
library(stringr)
library(plyr)

#Input File
Bedfile1='Path to/Data/Arm_bedGraph_fileList.txt.fivePoint_average.bed'
Bedfile2='Path to/Data/Norm_bedGraph_fileList.txt.fivePoint_average.bed'
AnnoFile='Path to/Data/Gene_Region.txt'
Outfile='Path to/Corss_Talk.pdf'





Bed=arrange(read.table(Bedfile1,sep='\t', col.names=c('Gene','X','Y'),colClasses = c("factor","numeric","numeric")),Gene)
Gene.x=aggregate(Bed$X,list(Bed$Gene),FUN=min)
rm(Bed)
Read_BedFile=function(fileName){
  Bed1=arrange(read.table(fileName,sep='\t', col.names=c('Gene','X','Y'),colClasses = c("factor","numeric","numeric")),Gene)
  Gene.num.Bed1=aggregate(Bed1$X,list(Bed1$Gene),FUN=length)
  Subtract=rep(Gene.x$x,Gene.num.Bed1$x)
  Bed1$X=Bed1$X-Subtract
  return(Bed1)
}

Bed1=Read_BedFile(Bedfile1)
Bed2=Read_BedFile(Bedfile2)


#Initialization
pdf(Outfile)
circos.par(start.degree=90,cell.padding = c(0.02, 0, 0.02, 0),gap.degree=c(rep(0.5,37),8))
circos.initialize(factors=Bed1$Gene, x = Bed1$X)

#First Cell
VSector=c('CD31_JAG1','CD4_FGF2','CD8_FGF2','DC_FGFBP3','Fib_JAG1','Mac_JAG1')

Colors=adjustcolor(c('darkorchid4','deeppink3','firebrick4','midnightblue','lightgoldenrod3','darkslategray4'),alpha.f=0.4)

circos.track(factors=Bed1$Gene,y=Bed1$Y,x=Bed1$X,bg.border="white",track.height=0.05,panel.fun=function(x,y){
  if (CELL_META$sector.index %in% VSector){
    circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2]+uy(5,'mm'),as.vector(str_split(CELL_META$sector.index,'_')[[1]])[1],cex=1)
  }
  circos.text(CELL_META$xcenter,CELL_META$ycenter,as.vector(str_split(CELL_META$sector.index,'_')[[1]])[2],cex=0.3)
})

Gene.type=unique(Bed1$Gene)
Cell.type=c()
for (gene.type in Gene.type){Cell.type=c(Cell.type,as.vector(str_split(gene.type,'_')[[1]])[1])}
GeneType=data.frame(GeneType=Gene.type,CellType=Cell.type)
Groups=aggregate(GeneType$GeneType,list(GeneType$CellType),length)
Gene.Groups=data.frame(row.names=unique(Cell.type))
ST=c(1)
EN=c(Groups$x[1])
for (i in Groups$x[-1]){
  ST=c(ST,EN[length(EN)]+1)
  EN=c(EN,EN[length(EN)]+i)
}

for (i in seq(length(ST))){
  draw.sector(get.cell.meta.data("cell.start.degree",sector.index=Gene.type[ST[i]]),
              get.cell.meta.data("cell.end.degree",sector.index=Gene.type[EN[i]]),
              rou1=get.cell.meta.data('cell.top.radius',track.index=1),
              rou2=get.cell.meta.data('cell.bottom.radius',track.index=1),
              col=Colors[i],
              border = Colors[i]
              )}
#Secode track(Arm ,red)
circos.track(factors=Bed1$Gene,y=Bed1$Y,x=Bed1$X,ylim=c(0,30),bg.border="black",bg.lwd=0.4,track.height=0.15,panel.fun=function(x,y){
  circos.lines(x,y,type='h',col='coral',border='coral')
  circos.axis(h='bottom',major.at=c(0,40000,80000,120000,160000,200000),labels=c('','','','',''),labels.cex=0.4,direction='inside')
})

circos.yaxis(sector.index='CD31_CCL21',track.index=2,labels.cex=0.5,labels=c('0','10','20','30'),at=c(0,10,20,30),side="left")


#Third track(Norm, darkblue)
circos.track(factors=Bed2$Gene,y=Bed2$Y,x=Bed2$X,ylim=c(0,30),bg.border="black",bg.lwd=0.4,track.height=0.15,panel.fun=function(x,y){
  circos.lines(x,y,type='h',col='royalblue4',border='royalblue4')
})
circos.yaxis(sector.index='CD31_CCL21',track.index=3,labels.cex=0.5,labels=c('0','10','20','30'),at=c(0,10,20,30),side="left")


#Add link
Region_table=arrange(read.table(AnnoFile,sep='\t',col.names=c('Gene','Start','End'),skip=1),Gene)
Region_table$Start=Region_table$Start-Gene.x$x
Region_table$End=Region_table$End-Gene.x$x

GetRegion=function(Sector){
  R=subset(Region_table,Gene==Sector)
  return(c(R$Start,R$End))
}

GenominLink=function(From,To,Color){
  for (i in seq(length(From))){
    sector1=From[i]
    sector2=To[i]
    circos.link(sector1,GetRegion(sector1),sector2,GetRegion(sector2),col=Color,border = NA)
  }
}

From1=c('DC_FAS','DC_FAS','DC_CCL22','DC_CCL22','DC_CD40','DC_FGFBP3','DC_FGFBP3','DC_FGFBP3','DC_FGFBP3','DC_FGFBP3','DC_ITGB6','DC_ITGB6','DC_ITGB6','DC_ITGB6')
To1=c('CD4_FASLG','CD8_FASLG','CD4_CCR4','CD8_CCR4','CD4_CD40LG','CD4_FGF2','CD8_FGF2','Fib_FGF2','Mac_FGF2','CD31_FGF2','CD4_TGFB1','CD8_TGFB1','Mac_TGFB1','CD31_TGFB1')
Color1=Colors[4]
From2=c('CD4_CCR7','CD4_CCR7','CD4_CCR7','CD4_DLL4','CD4_DLL4','CD4_DLL4','CD4_NOTCH1','CD4_NOTCH1','CD4_NOTCH1')
To2=c('Mac_CCL21','CD31_CCL21','Fib_CCL21','Mac_NOTCH4','CD31_NOTCH4','Fib_NOTCH4','Mac_JAG1','CD31_JAG1','Fib_JAG1')
Color2=Colors[2]
From3=c('CD8_CCR7','CD8_CCR7','CD8_CCR7','CD8_DLL4','CD8_DLL4','CD8_DLL4','CD8_NOTCH1','CD8_NOTCH1','CD8_NOTCH1')
To3=c('Mac_CCL21','CD31_CCL21','Fib_CCL21','Mac_NOTCH4','CD31_NOTCH4','Fib_NOTCH4','Mac_JAG1','CD31_JAG1','Fib_JAG1')
Color3=Colors[3]
From4=c('Mac_DLL4','Mac_DLL4')
To4=c('CD4_NOTCH4','CD8_NOTCH4')
Color4=Colors[6]
From5=c('CD31_DLL4','CD31_DLL4')
To5=c('CD4_NOTCH4','CD8_NOTCH4')
Color5=Colors[1]
GenominLink(From1,To1,Color1)
GenominLink(From2,To2,Color2)
GenominLink(From3,To3,Color3)
GenominLink(From4,To4,Color4)
GenominLink(From5,To5,Color5)

for (i in Gene.x$Group.1){
  circos.axis(sector.index = i,track.index=3,h='bottom',major.at=c(0,40000,80000,120000,160000,200000),labels=c('0b','40kb','80kb','120kb','',''),labels.cex=0.4,direction='inside')
}

dev.off()
