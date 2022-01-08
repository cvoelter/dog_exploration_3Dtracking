#function to calculates the effort per cell on a grid within which a
#transect (or any other straight line) was walked
#arguments are:
#grid.left: longitude of the left edge of the grid
#grid.right: longitude of the right edge of the grid
#grid.bottom: latitude of the bottom edge of the grid
#grid.top: latitude of the top edge of the grid
#nx: number of cells along the longitude
#ny: number of cells along the latitude
#xs: longitude of the transect's start
#xe: longitude of the transect's end
#ys: latitude of the transect's start
#ye: latitude of the transect's end
#all coordinates are assumed to be in UTM
#the function returns a list with two data frames
#$tr.segm comprises two columns showing the longitude (x) and latitude (y) of
#the points where the transect intersects with the grid lines as well as the
#coordinates of the start and end of the transect
#$effort comprises the effort per cell (in meters walked)
#cells are identified by their center points as well as by numbers consecutively
#labeling the cells along longitude and latitude from 1 to N
#last updated Oct. 8 2013
effort.per.cell<-function(grid.left, grid.right, grid.bottom, grid.top, nx, ny, xs, xe, ys, ye){
	slope=(ye-ys)/(xe-xs)
  intercept=ye-xe*slope
  x.vals=c()
  y.for.x=c()
  if(xs!=xe){
    x.vals=seq(grid.left, grid.right, length.out=nx+1)
    if(xs<=xe){
      x.vals=x.vals[x.vals>=xs & x.vals<=xe]
    }else{
      x.vals=x.vals[x.vals<=xs & x.vals>=xe]
    }
    y.for.x=intercept+slope*x.vals
  }
  y.vals=seq(grid.bottom, grid.top, length.out=ny+1)
  if(ys<=ye){
    y.vals=y.vals[y.vals>=ys & y.vals<=ye]
  }else{
    y.vals=y.vals[y.vals<=ys & y.vals>=ye]
  }
  if(xs!=xe){
    x.for.y=(y.vals-intercept)/slope
  }else{
    x.for.y=rep(xs, length(y.vals))
  }
  tr.segm=data.frame(x=c(xs, x.vals, x.for.y, xe), y=c(ys, y.for.x, y.vals, ye))
  tr.segm=tr.segm[order(tr.segm[,1], tr.segm[,2]),]
  x.width=(grid.right-grid.left)/nx
  y.width=(grid.top-grid.bottom)/ny
  x.name=1+((tr.segm$x[-1]+tr.segm$x[-nrow(tr.segm)])/2-grid.left)%/%x.width
  y.name=1+((tr.segm$y[-1]+tr.segm$y[-nrow(tr.segm)])/2-grid.bottom)%/%y.width
  x.name[x.name<1|x.name>nx]=NA
  y.name[y.name<1|y.name>ny]=NA
  effort=unlist(lapply(2:nrow(tr.segm), function(xxx){
    sqrt((tr.segm[xxx-1,1]-tr.segm[xxx,1])^2+(tr.segm[xxx-1,2]-tr.segm[xxx,2])^2)
  }))
  x.center=grid.left+x.width/2+(x.name-1)*x.width
  y.center=grid.bottom+y.width/2+(y.name-1)*y.width
  effort=data.frame(x.name, y.name, x.center, y.center, effort)
  return(list(tr.segm=tr.segm, effort=effort))
}
