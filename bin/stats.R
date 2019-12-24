data <- read.table("cputimes.csv",sep=";",header=T,dec=",")


output <- numeric()

for (cpu in c(1,2,4)) {
  for (alg in 1:2) {
    for (map in 1:4) {
      data.1 <- data[data$CPUs==cpu & data$ALG==alg & data$MAP==map,"TIME"]
   
      output <- c(output,cpu,alg,map,mean(data.1),sd(data.1))
    }
  }
}
output<- matrix(output,ncol=5,byrow=T)
colnames(output)<-c("CPUs","ALG","MAP","mean(TIME)","sd(TIME)")
output[["ALG"]]

output.alg2.map3 <- output[output[,"ALG"]==2 & output[,"MAP"]==3,]
output.alg2.map3
