
v.angle<-function(tailtip.loc, tailbase.loc, atlas.loc){
    tailtip.loc=tailtip.loc[c("x", "y")]
    tailbase.loc= tailbase.loc[c("x", "y")]
    atlas.loc=atlas.loc[c("x", "y")]
    tailtip.loc=tailtip.loc-tailbase.loc
    atlas.loc=atlas.loc-tailbase.loc
    result<-atan2(atlas.loc["y"],atlas.loc["x"]) - atan2(tailtip.loc["y"],tailtip.loc["x"]) 
    ind1 <- which(result > pi)
    ind2 <- which(result < -pi)
    result[ind1] <- result[ind1] - 2*pi
    result[ind2] <- result[ind2] + 2*pi
    return(result*180/pi) 
}

#tailtip.loc<-data.frame(x=1, y=0)
#atlas.loc<-data.frame(x=0, y=1)
#tailbase.loc<-data.frame(x=0, y=0)

#v.angle(tailtip.loc,tailbase.loc, atlas.loc)   

#left: negative values 
#right:positive values