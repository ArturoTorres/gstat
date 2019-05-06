## FNN local prediction
########################
library(sp)
library(spacetime)
library(gstat)
library(lattice)

# create n space-time points over [0,1] x [0,1] x [Now, Now+some days]
t0 = Sys.time() # now
n = 1e3
set.seed(13131) # fix outcomes
x = runif(n)
y = runif(n)
t = t0 + 1e6 * runif(n)
z = rnorm(n)
stidf = STIDF(SpatialPoints(cbind(x,y)), sort(t), data.frame(z=z))

stplot(stidf, number=21, main="random spatio-temporal noise")

# create a regular 20 x 20 x 20 grid of prediction locations:
grd = as(SpatialGrid(GridTopology(c(0.025,0.025), c(.05, .05), c(20,20))), "SpatialPixels")
tgrd = seq(min(t)+10000, max(t)-10000, length.out = 20)
stf = STF(grd, tgrd)

plot(index(stidf@time),1:n, main="sampled locations vs prediction steps")
abline(v=as.numeric(index(stf@time)), col="red")

# define a variogram model
sumMetricModel <- vgmST("sumMetric",
                        space=vgm(1/6, "Sph", 0.25, 1/60),
                        time =vgm(2/6, "Exp",  1e5, 1/60),
                        joint=vgm(0.4, "Exp", 0.3, 0.1),
                        stAni=1/1e6)
attr(sumMetricModel, "temporal unit") <- "secs"

dg <- data.frame(spacelag=rep(c(0.001,1:10)/10,6), 
                 timelag=rep(0:5*50e3, each=11))

wireframe(gamma ~ spacelag+timelag,
          variogramSurface(sumMetricModel, dist_grid = dg),
          scales=list(arrows=F),
          drape=T, col.regions=bpy.colors(),
          zlim=c(0,1.2),
          ylab="time [secs]",
          xlab="space [1]",
          main="imposed sum-metric model")

# space-time local kriging
time <- Sys.time()
locSpaceTimeKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmax=10)
time <- Sys.time() - time

stplot(locSpaceTimeKrig[,1:12], 
       main=paste("space-time local kriging", round(time,2), units(time)),
       col.regions=bpy.colors(), scales=list(draw=T))

# time slice wise local kriging
time <- Sys.time()
locTimeSliceKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmaxTime=c(-12,12)*3600)
time <- Sys.time() - time

stplot(locTimeSliceKrig[,1:12], 
       main=paste("time slice wise local kriging", round(time,2), units(time)),
       col.regions=bpy.colors(), scales=list(draw=T))

# time slice wise space-time local kriging
time <- Sys.time()
locTimeSliceSpaceTimeKrig <- krigeST(z~1, stidf, stf, sumMetricModel, nmaxTime=c(-12,12)*3600, nmax = 20)
time <- Sys.time() - time

stplot(locTimeSliceSpaceTimeKrig[,1:12], 
       main=paste("time slice wise space-time local kriging", round(time,2), units(time)),
       col.regions=bpy.colors(), scales=list(draw=T))
