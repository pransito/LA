### Grid package code selection

library(grid)

## defining the field of view (field of plotting)
## x and y define the center of the viewport
vp <- viewport(x=0.5,y=0.5,width=1, height=1)

## activation the new viewport

pushViewport(vp)

## coordinates lower left-hand corner (0,0); coordinates of upper right-hand corner (1,1)

pushViewport(vp)
# a rectangle (with dashed lines) on the border
# of the viewport:
grid.rect(gp=gpar(lty="dashed"))
# a circle centered at (.6,.4) with radius .3:
grid.circle(x=0.6, y=0.4, r=0.3)

## Consider the function stickperson() whose code is displayed
## here:
  
  stickperson <- function() {
    grid.circle(x=.5, y=.8, r=.1, gp=gpar(fill="yellow"))
    grid.lines(c(.5,.5), c(.7,.2)) # vertical line for body
    grid.lines(c(.5,.7), c(.6,.7)) # right arm
    grid.lines(c(.5,.3), c(.6,.7)) # left arm
    grid.lines(c(.5,.65), c(.2,0)) # right leg
    grid.lines(c(.5,.35), c(.2,0)) # left leg
  }

## The code uses the gp argument when drawing the circle;
## this argument takes a large number of graphical parameters
## which are specified by gpar(). gpar() is the grid version of R base par()

## multiple viewports

vp1 <- viewport(x=(1/3), y=1/3, width=0.6, height=0.3)
pushViewport(vp1)

## vp1 is now the child of vp because it got pushed into the viewport vp1
## if you want to go up to the parent viewport use upViewport()