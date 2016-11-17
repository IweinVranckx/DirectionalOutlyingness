
rm(list = ls())

source("DO_code.R")
load("DO_small_datasets.RData") # 1.3 Mb
load("DO_MRI_data.RData")       # 32 Mb
load("DO_Video_data.RData")     # 38 Mb

######################
# Family Income Data #
######################

temp=compScales(FamilyIncome)
sa=temp$sa
sb=temp$sb
med=temp$med
mad=mad(x = FamilyIncome)
sdovalues=SDO(FamilyIncome)$SDO
dovalues=DO_univ(FamilyIncome)[,3]
cutoffSDO=findoutliersSDO(x=FamilyIncome)$cutoff
cutoffDO=findoutliersDO(x=FamilyIncome)$cutoff

h=hist(FamilyIncome,breaks = seq(0,1e+7,by=5000),xlim=c(0,5e+5),
       plot=TRUE,cex.lab=1.5,ylab="",main="",
       xlab="Annual family income in USD")
title(ylab="Frequency", line=2.5, cex.lab=1.5)
arrows(x0 = med,y0 = 300,x1 = med-mad,y1=300,col="orange",length=0.1)
arrows(x0 = med,y0 = 300,x1 = med+mad,y1=300,col="orange",length=0.1)
points(x = med,y=300)
arrows(x0 = med,y0 = 200,x1 = med-sb,y1=200,col="blue",length=0.1)
arrows(x0 = med,y0 = 200,x1 = med+sa,y1=200,col="blue",length=0.1)
points(x = med,y=200)
legend("topright",legend = c("DO scales","SDO scales"),fill = 
         c("blue","orange"))
segments(x0 = (med+cutoffDO*sa),y0 = 0,x1 = (med+cutoffDO*sa),
         y1 = 100,col = "blue",lwd = 3)
segments(x0 = (med+cutoffSDO*mad),y0 = 0,x1 = (med+cutoffSDO*mad),
         y1 = 100,col = "orange",lwd = 3)


#################
# Bloodfat Data #
#################

axislengthx=seq(0,500,by=5)
axislengthy=seq(-200,1000,by=5)
xy=expand.grid(x=axislengthx,y=axislengthy)
xy=cbind(xy$x,xy$y)
nlev = seq(0.5,16,by=1)
z1=DO_multiv(X1 = Bloodfat,X2 = xy,algo="PP")$do
z1=matrix(data=z1,nrow=length(axislengthx),ncol=length(axislengthy))
z1[z1<1]=1
z2=SDO(x = Bloodfat,z = xy)$SDO
z2=matrix(data=z2,nrow=length(axislengthx),ncol=length(axislengthy))
z2[z2<1]=1

filled.contour(x=axislengthx,y=axislengthy,main="DO Contours",z=z1,
          color.palette =heat.colors,levels=nlev, plot.axes = 
          {points(Bloodfat[,1],Bloodfat[,2],pch=16);axis(1);axis(2)})

filled.contour(x=axislengthx,y=axislengthy,main="SDO Contours",z=z2,
          color.palette =heat.colors,levels=nlev, plot.axes = 
          {points(Bloodfat[,1],Bloodfat[,2],pch=16);axis(1);axis(2)})


##############
# Glass Data #
##############

weights=rep(1,750)
weights[1:13]=0

# first, treat glass data as univariate (not shown in paper):
ptm=proc.time()
DO_glass_uni=DO_functional(data=Glass[,,1],univ=TRUE,weights=weights)
print(proc.time()-ptm) # 1 second
FOM(DO_glass_uni,weights = weights)

DO_heatmap_sample(DO_glass_uni,weights = weights,XLAB="Wavelength",
                  breaks = c(0,50,100,150),labels = c(0,50,100,150))

# Now include the derivative, as in the paper:
ptm=proc.time()
DO_glass=DO_functional(data=Glass,weights=weights,algo="PP",
                       rmZeroes = TRUE)
print(proc.time()-ptm) # 100 seconds
FOM(DO_glass,weights = weights)

DO_heatmap_sample(DO_glass,weights = weights,XLAB="Wavelength",
                  breaks = c(0,50,100,150),labels = c(0,50,100,150))

############
# MRI Data #
############

# In the paper we analyzed the MRI data by the projection pursuit
# approach, which took a long time.
# To save computation time we use the componentwise procedure
# here. The resulting figures are almost identical.

ptm=proc.time()
DO_MRI_compwise=DO_functional(data=MRI,maxRatio = 2,algo="compWise")
print(proc.time()-ptm) # 50 seconds
FOM(DO_MRI_compwise)

DO_heatmap_image(t(apply(DO_MRI_compwise[387,,],2,rev)),cap=15)
DO_heatmap_image(t(apply(DO_MRI_compwise[92,,],2,rev)),cap=15)
DO_heatmap_image(t(apply(DO_MRI_compwise[126,,],2,rev)),cap=15)

############## 
# Video Data #
##############

ptm=proc.time()
DO_video_compwise=DO_functional(data=Video,algo="compWise",
                                rmZeroes = TRUE)
print(proc.time()-ptm) # 35 seconds
FOM(DO_video_compwise)

DO_heatmap_image(t(apply(DO_video_compwise[100,,],2,rev)),cap=75)
DO_heatmap_image(t(apply(DO_video_compwise[487,,],2,rev)),cap=75)
DO_heatmap_image(t(apply(DO_video_compwise[491,,],2,rev)),cap=75)
DO_heatmap_image(t(apply(DO_video_compwise[500,,],2,rev)),cap=75)

####################################################################

