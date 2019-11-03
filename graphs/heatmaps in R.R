## how to plot heatmaps in R
## see http://www.image.ucar.edu/GSP/Software/Fields/Help/image.plot.html

x<- 1:10; y<- 1:15; z<- outer( x,y,"+") 

image.plot(x,y,z) 
# or obj<- list( x=x,y=y,z=z); image.plot(obj)
# now add some points on diagonal with some clipping anticipated 
points( 5:12, 5:12, pch="X", cex=3)

# adding breaks and distinct colors
#        (lab.breaks are optional)

brk<- quantile( c(z))
image.plot(x,y,z, breaks=brk, col=rainbow(4), 
           lab.breaks=names(brk))

#
#fat (5 characters wide) and short (50% of figure)  color bar on the bottom
image.plot( x,y,z,legend.width=5, legend.shrink=.5, horizontal=TRUE) 

set.panel()

# Here is quick but quirky way to add a common legend to several plots. 
# The idea is leave some room in the margin and then over plot in this margin

par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel( 2,2) # 2X2 matrix of plots

# now draw all your plots using usual image command
for (  k in 1:4){
  image( matrix( rnorm(150), 10,15), zlim=c(-4,4), col=tim.colors())
}

par(oma=c( 0,0,0,1))# reset margin to be much smaller.
image.plot( legend.only=TRUE, zlim=c(-4,4)) 

# image.plot tricked into  plotting in margin of old setting 

set.panel() # reset plotting device