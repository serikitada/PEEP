setwd("~/Downloads/Sanger/01.Project/gitrepos/my_repositories/FlashFry/other_soft/DeepPE_modified/compare_with_web/")
web = read.table("web_example_DeepPE_Scores_woGtoC.txt",sep="\t",head=TRUE)
com = read.table("command_Scores.txt",sep="\t",head=TRUE)
com2 = read.table("command_Scores_woGtoC.txt",sep="\t",head=TRUE)
com = com[2:nrow(com),]
com2 =  com2[2:nrow(com2),]
View(web)
sum = cbind(web[,1:8],web[,10],web[,9],com)
sum2 =  cbind(web[,1:8],web[,10],web[,9],com2)
colnames(sum) = c("ID","Location","Wide","Guide","Edit","PBS_len","RTT_len","Extension","Model","Web","Command")
colnames(sum2) = c("ID","Location","Wide","Guide","Edit","PBS_len","RTT_len","Extension","Model","Web","Command")

View(sum)

cor(sum$Web,sum$Command)
dpe = sum[which(sum$Model=="DeepPE"),]
pos = sum[which(sum$Model=="PE_Position"),]
type = sum[which(sum$Model=="PE_type"),]
par(mfrow=c(2,1))
hist(dpe$Web,breaks=100,main="DeepPE Web score distribution",xlab="DeepPE Web score")
hist(dpe$Command,breaks=100,main="DeepPE Command score distribution",xlab="DeepPE Command score")
cor(dpe$Web,dpe$Command,method="spearman")
t.test(dpe$Web,dpe$Command,paired=TRUE,var.equal=FALSE)
plot(dpe$Web,dpe$Command)
par(new=T)
abline(a=0,b=1,lty=1,lwd=2)

par(mfrow=c(1,1))
p <-ggplot(dpe, aes(Web, Command)) + geom_point(size=0.8) +
  labs(title="DeepPE score comparison",x= "Web score", y = "Command score") + 
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14))
p + geom_abline(slope = 1,color="blue")

#t.test(dpe$Web, dpe$Command,paired=TRUE,var.equal=FALSE)
cor(pos$Web,pos$Command)
cor(type$Web,type$Command)
type$web
nrow(sum)
unique(sum[,3]) #25
table(sum[,3])


par(mfrow=c(1,1))
s = 1
spacer = unique(sum[,3])[s]
each = sum[which(sum[,3]==spacer),]
View(each)
colnames(each)
stripchart(Web ~ Model, data = each, 
           main = "Model",
           pch = 16,
           col = c("red", "green","blue"),
           group.names=c("DeepPE","PE_position","PE_type")
           vertical = TRUE)


#57
par(mfrow=c(5,5))
for (i in (1:5)){ #1:length(unique(sum[,3]))
  spacer = unique(sum[,3])[i]
  each = sum[which(sum[,3]==spacer),]
  stripchart(Web ~ Model, data = each, 
             main = "Model",
             ylab = "Web Score",
             pch = 16,
             col = c("red", "green","blue"),
             group.names=c("DeepPE","PE_position","PE_type"),
             vertical = FALSE)
}

sum10 = sum[which(sum[,3]%in%unique(sum[,3])[1:10]),]
p<-ggplot(sum10, aes(x=Wide, y=Web, color=Model, shape=Model)) +
  geom_jitter(position=position_dodge(0.8),size=2) +
  scale_color_brewer(palette="Blues",direction=-1) + theme_classic() +
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10"))+
  labs(y= "Web score", x = "Spacer")+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14))+
  theme(legend.title = element_text(size=14)) +
  theme(legend.text = element_text(size=14))
  
p

dpe$PBS_RTT_len = dpe$PBS_len + dpe$RTT_len
View(dpe)
ggplot(dpe, aes(x=Wide, y=Web, , size=2, color=PBS_RTT_len)) + 
  geom_point(size=1.5) + 
  scale_color_gradient(low="red", high="white")+
  scale_x_discrete(labels=c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","23","24","25"))+
  labs(y= "Web score", x = "Spacer")+
  theme(axis.title.y = element_text(size=14)) +
  theme(axis.title.x = element_text(size=14))+
  theme(legend.title = element_text(size=14)) +
  theme(legend.text = element_text(size=14))

out = c()
for (i in 1:25){
  out = c(out, paste('"',i,'"',sep=""))
}
out

# library
library(ggplot2)
# grouped boxplot
ggplot(sum, aes(x=Wide, y=Web, fill=Model)) + 
  geom_boxplot()



a <- runif(30,80,100)
b <-  runif(30,80,100)
c <-  runif(30,80,100)

plot(a, type="l", col="red1",ylim=c(70,100),lwd=2)
par(new=T)
plot(b, type="l", col="blue1",ylim=c(70,100),lwd=2)
par(new=T)
plot(c, type="l", col="green4",ylim=c(70,100),lwd=2)

d <- as.data.frame(rbind(a,b,c))
d[,31] = c("Set1","Set2","Set3")
colnames(d)[31] = "group"
d$group = c("Set1","Set2","Set3")
dd <- t(d)


?plot
ggplot( aes(x=, y=n, group=group, color=group)) +
  geom_line()