#Stem map should use smithsonian forestGeo data format
#Dead stems are removed because they do not contribute to richness
#Multiple stems are removed so they are only counted once. 
#Some ForestGeo data uses slightly different variable names, e.g. PX vs GX, or species vs latin. Some editing of the code may be necessary depending on the site data format. 

graphics.off(); rm(list=ls(all=TRUE)); library(vegan)
map <- read.table("U:\\Pudue Files\\Manuscripts\\Checkerboards\\Case Studies\\BCI\\Full.clean.BCI.map.txt",header=T,sep="\t",quote="")
map<-map[which(map$Status != "dead"),]; map<-map[which(map$Stem != "multiple"),]

i<-0; imax<-50

C.vec<-as.numeric(); C.sd.vec<-as.numeric(); C.vec.null<-as.numeric(); 
C.sd.null.vec<-as.numeric(); scale.vec.raw<-as.numeric()
min.plot.size<-1
max.plot.size<-sqrt(min(max(map$PX), max(map$PY))) 
plot.incriment<-1
scale<-min.plot.size

#####
#Triple loop iterates over scales, subsampling the map
######
while (scale<=max.plot.size)
	{
	S<-as.numeric()
	sp.list<-gsub(" ", ".", unique(map$species))
	comm.matrix<-rep(9999, length(sp.list))
	names(comm.matrix)<-sp.list
	comm.matrix<-data.frame(t(comm.matrix))[-1,]
	##########################
	#Loop to randomly sample areas
	# of size=scale (centimetres)
	#############################
	while (i<imax) 
		{
	
		x.lim<-seq(0, (max(map$PX)-scale), by=0.1)
		y.lim<-seq(0, (max(map$PY)-scale), by=0.1)

		x.min<-sample(x.lim, 1)
		y.min<-sample(y.lim, 1)
		x.max<-x.min+scale
		y.max<-y.min+scale

		map.min <- map[ which(map$PX > x.min & map$PY > y.min),]
		map.new <- map.min[ which(map.min$PX < x.max & map.min$PY < y.max),]
		
		sub.sp.list<-gsub(" ", ".", unique(map.new$species))
		S<-c(S, length(sub.sp.list))

			#########################
			#Loop to populate 
			#occurence matrix
			#########################
			j<-1
			k<-1
			temp.mat<-as.numeric() 
			while(j<=length(sub.sp.list))
				{ 
				m.n<-sub.sp.list[j]
				temp.vec<-as.numeric()
				while(k<=length(sp.list))
					{
					c.n<-sp.list[k]
					temp.vec<-c(temp.vec, (if(c.n==m.n) 1 else 0))
					k<-k+1
					}
				temp.mat<-rbind(temp.mat, temp.vec)
				k<-1	
				j<-j+1
				}
			if(length(sub.sp.list)==0) {x<-rep(0, length(sp.list))} else {x<-colSums(temp.mat)}
			comm.matrix<-rbind(comm.matrix, x)
			###########################
		j<-1
		i<-i+1
		}

	P<-length(sp.list)*(length(sp.list)-1)/2
	N<-imax*(imax-1)/2

	C.score<-sum(designdist(t(comm.matrix), method="((A-J)*(B-J))", terms="minimum"))/(N*P)	#Calculate mean Cij

######
#NULL MODEL averaging
######
	comm.matrix.null <- permatswap(comm.matrix, times = 1000, burnin = 30000, thin = 30000, mtype = "prab")

	 temp<-as.numeric
	for (i in 1:1000) 
		{ temp<-c(temp, (1/(N*P))*sum(designdist(t(comm.matrix.null$perm[[i]]), method="((A-J)*(B-J))", terms="minimum"))) }
			
	C.score.null<-mean(as.numeric(temp[-1]))
	C.sd.null<-sd(as.numeric(temp[-1]))

##########

	C.vec<-c(C.vec, C.score)
	C.vec.null<-c(C.vec.null, C.score.null)
	C.sd.null.vec<-c(C.sd.null.vec, C.sd.null)
	scale.vec.raw<-c(scale.vec.raw, scale)
	
	i<-0
	print(c(mean(C.vec), mean(C.vec.null), (scale^2)))
	flush.console()
	scale<-scale+plot.incriment
	names(comm.matrix)<-sp.list

	}

dev.new()
par(mfrow=c(2,1))
plot((scale.vec.raw^2)/10000, (C.vec), main="CASE STUDY NAME", col="red", 
	xlab="Quadrat size (Ha)", ylab="C-score")
points((scale.vec.raw^2)/10000, (C.vec.null))

effect.size<-(C.vec-C.vec.null)/(C.sd.null.vec)

plot((scale.vec.raw^2)/10000, effect.size, main="CASE STUDY NAME",  
	xlab="Quadrat size (Ha)", ylab="Effect size (z)")
lines(c(-1,5),c(2,2), col="red")


out<-cbind(scale.vec.raw^2/10000, C.vec, C.vec.null, C.sd.null.vec) 
#write.csv(out, "DIRECTORY") 
