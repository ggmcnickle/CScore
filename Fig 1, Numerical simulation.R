#Random community presence/absence matrix
#Figure 1
#McNickle and Lamb 
#Generates random area based samples for Checkerboard score
graphics.off(); rm(list=ls(all=TRUE))
library(ggplot2); library(vegan); library(gridExtra)

########################
########################
##  Parameters
########################
##  p = number of plots
##  s = number of species
##  a = area sampled (a from Kinsella was 0.1, 1 or 6 m^-2) 
##  B = Species area curve constant
##  z = species area curve exponent
########################
########################

x<-1000	

Num.C.Area<-function(p, B, z, a.max) 
{
	out<-as.numeric(); C.score<-as.numeric(); C.sd<-as.numeric()
	a<-as.numeric(); S<-as.numeric()
	for (i in 1:x)
		{
		a[i]<-a.max/(x/i); S[i]<-B*a[i]^z
		for (j in 1:p) 											#Loop 1 calculates Checkerboard scores
			{												#Loop2 generates a p*s matrix where p are samples, and s are species
			if (B*a[i]^z>=1) {temp.1 <- sample(1, B*a[i]^z, replace=TRUE)} else {temp.1 <- 0}					#Basic species area curve defines total number of presences
			if (length(temp.1<B*a.max^z)) {temp.0 <- sample(0, B*a.max^z-sum(temp.1), replace=TRUE)} else {temp.0 <- as.numeric()}		#non-presences are zeros
			temp   <- c(temp.1, temp.0)										#presences and absences combined	
			data <- sample(temp, B*a.max^z, replace=FALSE)									#presences and absences shuffled randomly
			out<-rbind(out, data)										#Bind new row to data set
			}
		C.score[i]<-mean(designdist(t(out), method="((A-J)*(B-J))"))/(p*(p-1)/2)	#Calculate mean Cij
		C.sd[i]<-sd(designdist(t(out), method="((A-J)*(B-J))"))	#StDev of Cij
		out<-as.numeric()										#reset for loop i
		}
out<-cbind(a,S, C.score, C.sd)
out
}

z<-0.1; z1<-cbind(rep(z,x), Num.C.Area(50, 25, z, 100)) 
z<-0.2; z2<-cbind(rep(z,x), Num.C.Area(50, 25, z, 100)) 
z<-0.3; z3<-cbind(rep(z,x), Num.C.Area(50, 25, z, 100)) 
z<-0.4; z4<-cbind(rep(z,x), Num.C.Area(50, 25, z, 100)) 
z.panel<-data.frame(rbind(z1, z2, z3, z4)); colnames(z.panel)[1]<-"z"

B<-10; B1<-cbind(rep(B,x), Num.C.Area(50, B, 0.25, 100)) 
B<-20; B2<-cbind(rep(B,x), Num.C.Area(50, B, 0.25, 100)) 
B<-30; B3<-cbind(rep(B,x), Num.C.Area(50, B, 0.25, 100)) 
B<-40; B4<-cbind(rep(B,x), Num.C.Area(50, B, 0.25, 100)) 
B.panel<-data.frame(rbind(B1, B2, B3, B4)); colnames(B.panel)[1]<-"B"

panel.1a<-ggplot(z.panel, aes(a, C.score, colour=as.factor(z)) ) + geom_line() +
	labs(x="Sample unit scale", y="C-score") + theme_bw() + scale_colour_discrete(name="  z") 
panel.1b<-ggplot(z.panel, aes(a, S, colour=as.factor(z)) ) + geom_line() +
	labs(x="Sample unit scale", y="Species richness") + theme_bw() + scale_colour_discrete(name="  z") 
panel.1c<-ggplot(B.panel, aes(a, C.score, colour=as.factor(B)) ) + geom_line() +
	labs(x="Sample unit scale", y="C-score") + theme_bw() + scale_colour_discrete(name="  B") 
panel.1d<-ggplot(B.panel, aes(a, S, colour=as.factor(B)) ) + geom_line() +
	labs(x="Sample unit scale", y="Species richness") + theme_bw() + scale_colour_discrete(name="  B") 

grid.arrange(panel.1a, panel.1b, panel.1c, panel.1d, ncol=2)

 
