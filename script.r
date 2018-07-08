############################
# Projet credit
############################

# Library

# install.packages("latex2exp")
library(latex2exp)
library(ggplot2)

# Import clean data

data.raw <- read.csv("data.csv",header=TRUE)

entire.Date <- factor(data.raw[,1],levels=c(levels(data.raw[,1]),""))
data.final <- data.frame(Date=entire.Date)

for (i in seq(1,length(data.raw),2)) {
	temp <- factor(data.raw[,i],levels=levels(entire.Date))
	idx <- which(entire.Date %in% temp)
	data.temp <- rep(NA,length(entire.Date))
	data.temp[idx] <- data.raw[1:length(idx),i+1]

	data.final <- data.frame(data.final,data.temp)
}

data.final$Date <- as.Date(data.final$Date,format="%m/%d/%Y")

names <- names(data.raw)[seq(1,length(data.raw),2)]
names <- c("Date",names)
names(data.final) <- names

data <- data.final
rm(data.raw,data.final,temp,idx,data.temp,entire.Date)


# Parameters estimation

end <- length(data$Date)
period <- 5

yield <- data.frame(data$Date[(period+1):end])
for (i in 2:length(data)) {
	temp <- rep(NA,end-period)	
	for (j in 1:(end-period)) {
		temp[j] <- data[j+period,i] / data[j,i]
	}
	yield <- data.frame(yield,temp)
}
names(yield) <- names

mu <- colMeans(yield[,2:length(yield)],na.rm=TRUE)
omega <- cov(yield[,2:length(yield)],use="pairwise.complete.obs")

r <- -40*1e-4
v <- 1


# Verify no only NA

fun <- function(x) sum(is.na(x))
temp <- sapply(yield,fun)==length(yield$Date)
sum(temp)


# Efficient frontier (with NRA)

one <- rep(1,length(mu))
inv.omega <- solve(omega)

A <- as.numeric(mu %*% inv.omega %*% mu)
B <- as.numeric(mu %*% inv.omega %*% one)
C <- as.numeric(one %*% inv.omega %*% one)


std <- seq(0,20*1e-4,1e-7)
ef.nra <- v*(1+r) + std*sqrt(A+(1+r)**2*C-2*(1+r)*B)

lambda.ra.plus <- sqrt((A*C-B**2)/(std**2*C-v**2))
phi.ra.plus <- (B-lambda.ra.plus*v)/C
ef.ra.plus <- (A-phi.ra.plus*B)/lambda.ra.plus

lambda.ra.minus <- -sqrt((A*C-B**2)/(std**2*C-v**2))
phi.ra.minus <- (B-lambda.ra.minus*v)/C
ef.ra.minus <- (A-phi.ra.minus*B)/lambda.ra.minus

std.min <- v/sqrt(C)
ef.ra.min <- v*B/C

idx.diff <- which.min(abs(ef.nra-ef.ra.plus))
std.diff <- std[idx.diff]
ef.diff <- ef.nra[idx.diff]

xdist <- 2
ydist <- 5e-4

pdf(file="./figure/ef.pdf")
plot(std*1e4,ef.nra,
	type="l",
	lwd=2,
	xlab=TeX("$\\sqrt{Var\\[V_1\\]}$ (bp)"),
	ylab=TeX("$E\\[V_1\\]$ (ccy)"))
lines(std*1e4,ef.ra.plus,
	col="blue",
	lwd=2)
lines(std*1e4,ef.ra.minus,
	col="blue",
	lwd=2,
	lty=2)
rect(std.min*1e4-xdist,ef.ra.min-ydist,std.min*1e4+xdist,ef.ra.min+ydist,
	lty=2)
grid()
dev.off()


# Zoom in

pdf(file="./figure/ef_zoom.pdf")
plot(std*1e4,ef.nra,
	type="l",
	lwd=2,
	xlim=c(std.min*1e4-xdist,std.min*1e4+xdist),
	ylim=c(ef.ra.min-ydist,ef.ra.min+ydist),
	xlab=TeX("$\\sqrt{Var\\[V_1\\]}$ (bp)"),
	ylab=TeX("$E\\[V_1\\]$ (ccy)"))
lines(std*1e4,ef.ra.plus,
	col="blue",
	lwd=2)
lines(std*1e4,ef.ra.minus,
	col="blue",
	lwd=2,
	lty=2)
points(std.min*1e4,ef.ra.min,
	col="green",
	pch=3,
	lwd=2)
points(std.diff*1e4,ef.diff,
	col="red",
	pch=3,
	lwd=2)
legend("topleft",inset=0.02,
	c("With NRA","Whithout NRA","Min variance portfolio","Same portfolio"),
	col=c("black","blue","green","red"),
	lwd=rep(2,4),
	lty=c(1,1,NA,NA),
	pch=c(NA,NA,3,3))
grid()
dev.off()



# Estimation matrix

var.diag <- sqrt(diag(omega))
correlation <- omega / (var.diag%*%t(var.diag))
# corr <- cor(yield[,2:length(yield)],use="pairwise.complete.obs")

eigen.values <- eigen(correlation)[[1]]

Q <- dim(yield)[1]/(dim(yield)[2]-1)
ev <- seq(0, 30, 0.05)

# lambda.max <- 1 + 1/Q + 2*sqrt(1/Q)
# lambda.min <- 1 + 1/Q - 2*sqrt(1/Q)

# rho <- ifelse(ev>lambda.min & ev<lambda.max, Q/2/pi*sqrt((lambda.max-ev)*(ev-lambda.min))/ev, 0)

N <- length(eigen.values)
sigma <- 1
rho <- matrix(data=NA, nrow=length(ev), ncol=N+1)

y <- sigma^2
lambda.max <- y*(1 + 1/Q + 2*sqrt(1/Q))
lambda.min <- y*(1 + 1/Q - 2*sqrt(1/Q))
rho[,1] <- ifelse(ev>lambda.min & ev<lambda.max, Q/2/pi/y*sqrt((lambda.max-ev)*(ev-lambda.min))/ev, 0)
for (i in 1:N) {
	y <- sigma^2 - sum(eigen.values[1:i])/N
	lambda.max <- y*(1 + 1/Q + 2*sqrt(1/Q))
	lambda.min <- y*(1 + 1/Q - 2*sqrt(1/Q))
	rho[,i+1] <- ifelse(ev>lambda.min & ev<lambda.max, Q/2/pi/y*sqrt((lambda.max-ev)*(ev-lambda.min))/ev, 0)
}

# Hist eigenvalues

# pdf(file="./figure/hist_ev.pdf")
a <- hist(eigen.values,
	breaks=50,
	freq=TRUE,
	main=NULL,
	xlab="Valeurs propres",
	ylab="Fréquence",
	col="darkblue")
grid()
# dev.off()


# Add dist M-P

h <- hist(eigen.values[eigen.values<2],
	breaks=100,
	probability=TRUE,
	main=NULL,
	xlab="Valeurs propres",
	ylab="Probabilité",
	col="darkblue")
grid()
lines(ev,rho[,1],
	lwd=2,
	col="red")
lines(ev,rho[,2],
	lwd=2,
	col="green")
lines(ev,rho[,4],
	lwd=2,
	col="violet")

gg <- ggplot(data.frame(eigen.values),aes(eigen.values)) +
	geom_density(adjust=0.5)
gg_b <- ggplot_build(gg)
plot(gg_b$data[[1]]$x,gg_b$data[[1]]$y,
	type="l")
grid()
lines(ev,rho[,1],
	lwd=2,
	col="red")
lines(ev,rho[,2],
	lwd=2,
	col="green")
lines(ev,rho[,4],
	lwd=2,
	col="violet")




# Test

x1 <- 1:5
y1 <- 2*x1
y2 <- 3*x1
t <- gather(data.frame(x1,y1,y2),variable,value,-x1)
ggplot(t,aes(x1,value,color=variable)) + geom_line()+scale_colour_manual(values=c("black", "orange"))


