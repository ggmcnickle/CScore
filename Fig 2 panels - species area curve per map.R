#Stem map should use smithsonian forestGeo data format
#Dead stems are removed because they do not contribute to richness
#Multiple stems are removed so they are only counted once. 
#Some ForestGeo data uses slightly different variable names, e.g. PX vs GX, or species vs latin. Some editing of the code may be necessary depending on the site data format. 

graphics.off(); rm(list=ls(all=TRUE))
library(ggplot2); library(reshape); library(gridExtra)

map <- read.table("U:\\Pudue Files\\Manuscripts\\Checkerboards\\Case Studies\\BCI\\Full.clean.BCI.map.txt",header=T,sep="\t",quote="")
map<-map[which(map$Status != "dead"),]; map<-map[which(map$Stem != "multiple"),]
scale<-1; increment<-1; effort<-50
max.plot.size<-min(max(map$PX), max(map$PY))

sp.list<-unique(map$species)	
	comm.matrix<-rep(9999, length(sp.list))
	names(comm.matrix)<-sp.list
	comm.matrix<-data.frame(t(comm.matrix))[-1,]
	i<-0
	imax<-effort
	S<-as.numeric()
	scale.vec<-as.numeric()
	max.quadrat<-min(c(max(map$PX)+increment,max(map$PY)+increment))
	while (scale <= max.plot.size^2)
		{
		while (i<imax) 
			{
			
			x.lim<-seq(0, (max(map$PX))-scale, by=increment)
			y.lim<-seq(0, (max(map$PY))-scale, by=increment)

			x.min<-sample(x.lim, 1)
			y.min<-sample(y.lim, 1)
			x.max<-x.min+scale
			y.max<-y.min+scale

			map.min <- map[ which(map$PX > x.min & map$PY > y.min),]
			map.new <- map.min[ which(map.min$PX < x.max & map.min$PY < y.max),]
		
			sub.sp.list<-gsub(" ", ".", unique(map.new$species))
			S<-c(S, length(sub.sp.list))
			scale.vec<-c(scale.vec,(scale^2)/10000)

			i<-i+1
			}
		i<-0
		scale<-scale+increment
		}

result<-data.frame(cbind(scale.vec,S))
#write.csv(result, "FILE LOCATION AND NAME", row.names=F) #insert file location and name to save output

panel.a<-ggplot(result, aes(scale.vec, S) ) + geom_point(alpha=.1) +
	labs(x="Quadrat size (Ha)", y="Cumulative species richness") + geom_smooth( method = 'loess') +
	theme_bw() 

	
plot(panel.a)

SAR.nls <- nls(S ~ a*(scale.vec*.0001)^b, data=result,
 start=list('a'=30,  'b'=.1))

SAR.nls
