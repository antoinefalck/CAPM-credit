###################################
############################
#
#		 Projet credit
#
############################
###################################



###################################
# Library

# install.packages("latex2exp")
# install.packages("tidyverse")
library(latex2exp)
library(ggplot2)


###################################
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



###################################
# Estimation matrix

var.diag <- sqrt(diag(omega))
correlation <- omega / (var.diag%*%t(var.diag))

eigen.values <- eigen(correlation)[[1]]
eigen.vectors <- eigen(correlation)[[2]]

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

# Or with ggplot2
# ggplot(data.frame(eigen.values),aes(x=eigen.values,y=..density..)) +
# 	geom_histogram(fill="darkblue")


# MP distribution and find lambda star

gg <- ggplot(data.frame(eigen.values),aes(eigen.values)) +
	geom_density(adjust=0.3)
gg_b <- ggplot_build(gg)
gg.x <- gg_b$data[[1]]$x
gg.y <- gg_b$data[[1]]$y


Q <- dim(yield)[1]/(dim(yield)[2]-1)
N <- length(eigen.values)
sigma <- 1
rho <- matrix(data=NA, nrow=length(gg.x), ncol=N+1)
mse <- rep(NA,N+1)

y <- sigma^2
lambda.max <- y*(1 + 1/Q + 2*sqrt(1/Q))
lambda.min <- y*(1 + 1/Q - 2*sqrt(1/Q))
idx <- gg.x>lambda.min & gg.x<lambda.max
rho[,1] <- ifelse(idx, Q/2/pi/y*sqrt((lambda.max-gg.x)*(gg.x-lambda.min))/gg.x, 0)
mse[1] <- mean((rho[idx,1]-gg.y[idx])^2)

for (i in 1:N) {
	y <- sigma^2 - sum(eigen.values[1:i])/N
	lambda.max <- y*(1 + 1/Q + 2*sqrt(1/Q))
	lambda.min <- y*(1 + 1/Q - 2*sqrt(1/Q))
	idx <- gg.x>lambda.min & gg.x<lambda.max

	rho[,i+1] <- ifelse(idx, Q/2/pi/y*sqrt((lambda.max-gg.x)*(gg.x-lambda.min))/gg.x, 0)
	mse[i+1] <- mean((rho[idx,i+1]-gg.y[idx])^2)
}


# Plot smoothed density and MP distribution

# pdf(file="./figure/smoothed_hist_ev_and_mp.pdf")
plot(gg.x[gg.x<2],gg.y[gg.x<2],
	type="l",
	ylim=c(0,3.3),
	lwd=2,
	xlab="Valeurs propres",
	ylab="Densité")
grid()
lines(gg.x,rho[,1],
	lwd=2,
	col="darkred")
lines(gg.x,rho[,2],
	lwd=2,
	col="darkgreen")
lines(gg.x,rho[,4],
	lwd=2,
	col="darkviolet")
legend("topleft",inset=0.02,
	c("Densité lissée",
		paste("M-P dist. avec y =",round(sigma^2,2)),
		paste("M-P dist. avec y =",round(sigma^2 - sum(eigen.values[1:1])/N,2)),
		paste("M-P dist. avec y =",round(sigma^2 - sum(eigen.values[1:3])/N,2))),
	col=c("black","darkred","darkgreen","darkviolet"),
	lwd=rep(2,4))
# dev.off()


# MSE like in algo 1

# pdf(file="./figure/eqm.pdf")
plot(0:6,mse[1:7],
	xlab="Index",
	ylab="Erreur quadratique moyenne",
	lwd=2)
grid()
# dev.off()


# New correlation matrix

Q <- dim(yield)[1]/(dim(yield)[2]-1)
N <- length(eigen.values)
sigma <- 1
lambda.max <- sigma^2*(1 + 1/Q + 2*sqrt(1/Q))
# lambda.max <- 3*lambda.max

L <- sum(eigen.values<lambda.max)
ev.const <- (N*sigma^2-sum(eigen.values[eigen.values>=lambda.max])) / L
eigen.values.new <- rep(ev.const,N)
eigen.values.new[1:(N-L)] <- eigen.values[1:(N-L)]

correlation.new <- eigen.vectors %*% diag(eigen.values.new) %*% solve(eigen.vectors)
diag(correlation.new) <- rep(1,N)
omega.new <- correlation.new * (var.diag%*%t(var.diag))


# Efficient frontiers new omega

inv.omega.new <- solve(omega.new)

A.new <- as.numeric(mu %*% inv.omega.new %*% mu)
B.new <- as.numeric(mu %*% inv.omega.new %*% one)
C.new <- as.numeric(one %*% inv.omega.new %*% one)

ef.nra.new <- v*(1+r) + std*sqrt(A.new+(1+r)**2*C.new-2*(1+r)*B.new)

lambda.ra.plus.new <- sqrt((A.new*C.new-B.new**2)/(std**2*C.new-v**2))
phi.ra.plus.new <- (B.new-lambda.ra.plus.new*v)/C.new
ef.ra.plus.new <- (A.new-phi.ra.plus.new*B.new)/lambda.ra.plus.new

# lambda.ra.minus <- -sqrt((A*C-B**2)/(std**2*C-v**2))
# phi.ra.minus <- (B-lambda.ra.minus*v)/C
# ef.ra.minus <- (A-phi.ra.minus*B)/lambda.ra.minus

# std.min <- v/sqrt(C)
# ef.ra.min <- v*B/C

# idx.diff <- which.min(abs(ef.nra-ef.ra.plus))
# std.diff <- std[idx.diff]
# ef.diff <- ef.nra[idx.diff]

# xdist <- 2
# ydist <- 5e-4

plot(std*1e4,ef.nra,
	type="l",
	lwd=2)
lines(std*1e4,ef.nra.new,
	col="blue",
	lwd=2)
lines(std*1e4,ef.ra.plus)
lines(std*1e4,ef.ra.plus.new,
	col="blue")

