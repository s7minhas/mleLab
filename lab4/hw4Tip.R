# Clear workspace
rm(list=ls())

# Load map
load('~/Downloads/crisp.map.rda')

# Set up map for plotting
## Lets also adjust the colors for the country borders in the map, this is 
## easily done by setting the border parameter to grey
plot(crisp.map, border='grey')

# Peak at one of the data structures included in the spatial polygon
attributes(crisp.map)$data[1:10,]

# pull out capital city lat and long
capLongLat=attributes(crisp.map)$data[,c('CNTRY_NAME', 'CAPNAME', 'CAPLONG', 'CAPLAT', 'COWCODE')]
capLongLat[1:10,]

# Now lets say that we want to add points to represent the locations of
# every capital city in the icews dataframe
## All we would need to do is reference the longitude (x position) and 
## latitude (y position) of every point. 
points(capLongLat$CAPLONG, capLongLat$CAPLAT)

### What this should tell you is that this map we've plotted is in
### many ways no different than any simple R plot. The only real 
### difference so far is that the x axis now represents longitude
### and the y latitude. 

# Now lets say we want to show that Russia "sent" a conflict to China
# First lets just pull out the relevant data points
chinaRussia=capLongLat[which(capLongLat$CNTRY_NAME %in% c('China', 'Russia')),]
chinaRussia

# Now we are going to use the arrows function to create a line that shows 
# Russia sent a conflict to China. 
## To do this we also have to choose a point for which to start and end the line, 
## lets just choose their capitals.
### To show you guys the exact inputs for this arrows function, I'm going to create
### a few more objects:
russLon=chinaRussia[1,3]
russLat=chinaRussia[1,4]
chinLon=chinaRussia[2,3]
chinLat=chinaRussia[2,4]

# Additionally, I'm going to set a few extra parameters. length controls
# the size of the arrow and lwd controls the line width.
arrows(x0=russLon, y0=russLat, x1=chinLon, y1=chinLat, length=.1, lwd=2)

# The points are just getting in the way so lets replot without them
plot(crisp.map, border='grey')
arrows(x0=russLon, y0=russLat, x1=chinLon, y1=chinLat, length=.1, lwd=2)