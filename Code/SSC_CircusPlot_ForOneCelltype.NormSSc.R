Args <- commandArgs()
print(Args)

library(circlize)
library(stringr)
library(plyr)


inFile=Args[6]


Remove0=function(V){
  newV=c()
  for (i in V){
    if (i!='0' & i!='0.0'){newV=c(newV,i)}
  }
  return(newV)
}


Remove1=function(chr){
  if (str_sub(chr,1,1)=='+'){return(c('+',str_sub(chr,2,)))}
  if (str_sub(chr,1,1)!='+'){return(c('-',str_sub(chr,2,)))}
}

RemoveJ=function(V){
  VN=c()
  for (chr in V){
    VN=c(VN,Remove1(chr)[2])
  }
  return(VN)
}


################################## Adjust DataFrame to fit the circosplot
#as.character(ColorDF['CD4',])
#adjustcolor(c('darkslategray4'),alpha.f=0.4)
ColorDF=data.frame(color=c('#EA2929','#F47D32','#6FBE47','#3D55A4'),bgcolor=adjustcolor(c('#EA2929','#F47D32','#6FBE47','#3D55A4'),alpha.f=0.4),linecolor=adjustcolor(c('#EA2929','#F47D32','#6FBE47','#3D55A4'),alpha.f=0.8),row.names = c('CD4','CD8','DC','Fib'))

NX=read.table(inFile,header = TRUE,sep='\t')
AllCells=colnames(NX)
OtherCells=AllCells[seq(2,length(AllCells))]
CellType=AllCells[1]

DCGenes=sort(unique(Remove0(NX$DC)))
CD4Genes=sort(unique(Remove0(NX$CD4)))
CD8Genes=sort(unique(Remove0(NX$CD8)))
FibGenes=sort(unique(Remove0(NX$Fib)))

Factors=as.factor(c(rep('DC',length(DCGenes)),rep('CD4',length(CD4Genes)),rep('CD8',length(CD8Genes)),rep('Fib',length(FibGenes))))


XL=c(seq(1,length(DCGenes)),seq(1,length(CD4Genes)),seq(1,length(CD8Genes)),seq(1,length(FibGenes)))
YL=c(rep(1,length(DCGenes)),rep(1,length(CD4Genes)),rep(1,length(CD8Genes)),rep(1,length(FibGenes)))


GenesL=c(DCGenes,CD4Genes,CD8Genes,FibGenes)

Circos.DF=data.frame(factor=Factors,X=XL,Y=YL,Genes=GenesL)

XlimM=matrix(c(0,0,0,0,length(DCGenes)+1,length(CD4Genes)+1,length(CD8Genes)+1,length(FibGenes)+1),ncol=2)
row.names(XlimM)=c('DC','CD4','CD8','Fib')

################################# initialize Circos
Outfile=paste0(inFile,'.CircosPlot.pdf')
pdf(Outfile,width=5.5,height=5.5)
circos.par(start.degree=90,cell.padding = c(0, 3, 0, 3),gap.degree=3)
circos.initialize(factors=Circos.DF$factor, x = Circos.DF$X,xlim=XlimM)

#############################first circle
circos.track(factors=Circos.DF$factor,x=Circos.DF$X,y=Circos.DF$Y,ylim=c(0,1),bg.border="white",track.height=0.2,panel.fun=function(x,y){
  circos.text(CELL_META$xcenter,CELL_META$cell.ylim[2]+uy(5,'mm'),CELL_META$sector.index,cex=1.0,col=as.character(ColorDF[CELL_META$sector.index,'color']))
  #Xlim=as.integer(CELL_META$cell.xlim)
  Xlim=XlimM[CELL_META$sector.index,]
  Sector.index=get.cell.meta.data("sector.index")
  Labels=RemoveJ(as.vector(Circos.DF[Circos.DF$factor==Sector.index,]$Genes))
  print(Xlim)
  
  draw.sector(get.cell.meta.data("cell.start.degree",sector.index=CELL_META$sector.index),
              get.cell.meta.data("cell.end.degree",sector.index=CELL_META$sector.index),
              rou1=get.cell.meta.data('cell.top.radius',track.index=1),
              rou2=get.cell.meta.data('cell.bottom.radius',track.index=1),
              col=as.character(ColorDF[CELL_META$sector.index,'bgcolor']),
              border =as.character(ColorDF[CELL_META$sector.index,'bgcolor']))
             
  circos.axis(h='bottom',major.at=seq(1,Xlim[2]-0.5),labels=Labels,labels.facing='reverse.clockwise',direction='outside',minor.ticks=0,
              labels.cex=0.5)
})

########################## Add Link
Rows=row.names(NX)
DrawLink=function(S1,S2){
  Sf=NX[S1]
  St=NX[S2]
  for (i in Rows){
      Fg=as.character(Sf[i,])
      Tg=as.character(St[i,])
      fi=Circos.DF[(Circos.DF$factor==S1)&(Circos.DF$Genes==Fg),]$X
      ti=Circos.DF[(Circos.DF$factor==S2)&(Circos.DF$Genes==Tg),]$X
      if ((Tg!='0') & (Tg!='0.0')){
        Fg=Remove1(Fg)
        Tg=Remove1(Tg)
        if ((Fg[1]=='+')&(Tg[1]=='+')){
          print('+,+')
          circos.link(S1,fi,S2,ti,h=1,col='firebrick',lwd=1,arr.length=0.2,directional=0)
        }else if((Fg[1]=='-')&(Tg[1]=='-')){
          circos.link(S1,fi,S2,ti,h=1,col='steelblue3',lwd=1,directional=0,arr.length=0.2)
        
        }
      }
  }
}

for (ct in OtherCells){DrawLink(CellType,ct)}

dev.off()



  
  
 
  









