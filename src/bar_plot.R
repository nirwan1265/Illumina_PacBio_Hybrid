# Data frame
data=data.frame(ACE1_dominance=3.211,ACE1_epistasis=6.678,ACE_0.2X43_X_N1111=-0.0074,ACE_0.2X43_X_N2833=-0.00389,ACE_0.2X43_X_N1009=0.00328,ACE_0.2X43_X_N1124=-0.0102,ACE_0.2X43_X_N2543=0.0169)
#data=data.frame(ACE1_dominance=c(3.211,"dominant"),ACE1_epistasis=c(6.678,"recessive"),ACE_0.2X43_X_N1111=c(-0.0074,"recessive"),ACE_0.2X43_X_N2833=c(-0.00389,"recessive"),ACE_0.2X43_X_N1009=c(0.00328,"recessive"),ACE_0.2X43_X_N1124=c(-0.0102,"recessive"),ACE_0.2X43_X_N2543=c(0.0169,)
# Group information
#groups <- c("dominant", "recessive", "recessive", "recessive", "recessive", "recessive", "recessive")

# Creating boundaries
lower=c(-0.015,0.015)
upper=c(1,8)
y_outer=21

# Creating spans
lowspan=c(0,11)
topspan=c(lowspan[2]+1,21)


# Figure labels
ylabel="y-axis value"
xlabel="x-axis value"
legendtext=c("ACE1_dominance","ACE1_epistasis","ACE_0.2X43_X_N1111","ACE_0.2X43_X_N2833","ACE_0.2X43_X_N1009","ACE_0.2X43_X_N1124","ACE_0.2X43_X_N2543")

# Functions
cnvrt.coords <-function(x,y=NULL){
  xy <- xy.coords(x,y, recycle=TRUE)
  cusr <- par('usr')
  cplt <- par('plt')  
  plt <- list()
  plt$x <- (xy$x-cusr[1])/(cusr[2]-cusr[1])
  plt$y <- (xy$y-cusr[3])/(cusr[4]-cusr[3])
  fig <- list()
  fig$x <- plt$x*(cplt[2]-cplt[1])+cplt[1]
  fig$y <- plt$y*(cplt[4]-cplt[3])+cplt[3]
  return( list(fig=fig) )
}

subplot <- function(fun, x, y=NULL){
  old.par <- par(no.readonly=TRUE)
  on.exit(par(old.par))
  xy <- xy.coords(x,y)
  xy <- cnvrt.coords(xy)$fig
  par(plt=c(xy$x,xy$y), new=TRUE)
  fun
  tmp.par <- par(no.readonly=TRUE)
  return(invisible(tmp.par))
}

# Will print in another window
quartz()

# Plotting
plot(c(-0.5,1),c(0,y_outer),type='n',axes=FALSE,ylab=ylabel,xlab='',lwd=7)


# Color palette for the bars
#colors <- c("#FF7F50", "#87CEFA")

subplot({
  y <- as.matrix(data)
  bp <- barplot(y,col=heat.colors(2),ylim=lower,xpd=FALSE,las=3)
},x=c(0,1),y=lowspan)

subplot({
  bp <- barplot(y, col=heat.colors(2), ylim=upper, xpd=FALSE,
                names.arg=vector(mode="character",length=length(data)))
}, x=c(0,1), y=topspan)


