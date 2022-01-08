dog.view.angle.and.eye.pos<-function(head.pos, snout.pos, eye.frac, eye.angle){
	##determines the viewing direction from an eye from the position of the top of the head and the snout,
		##the fraction of the distance between top of head and snout where the eye is located (0=eye at top of head,
			##1=eye at snout), and the vertical angle between the line connecting the top of the head and the snout and the
			##viewing direction
	##also determines the position of the eye
	##head.pos and snout.pos indicate the position of the top of the head and the snout in a cartesian coordinate system
		##must be named vectors; names must be x, y, and z
		##x and y denote horizontal (E-W) and vertical (N-S) position, and z the elevation
	##eye.frac scalar from interval [0, 1]
	##eye.angle vertical angle (in radians) between the line connecting the top of the head and the snout and the
		##viewing direction (must be positive)
	##returns the looking direction (in radians) with h being the horizontal (x-y) and v being the vertical (z) looking direction
	##the eye position returned is indicated in x, y, z-coordinates
	head.pos=head.pos[c("x", "y", "z")]
	snout.pos=snout.pos[c("x", "y", "z")]
	snout.pos=snout.pos-head.pos
	angle=c(h=as.vector(atan(snout.pos["y"]/snout.pos["x"])), v=NA)
	angle["v"]=as.vector(atan(snout.pos["z"]/sqrt(sum(snout.pos[c("x", "y")]^2))))
	#browser()
	angle["v"]=angle["v"]+eye.angle
	if(angle["v"]>pi/2){
		angle["v"]=pi-angle["v"]
	}
	if(angle["v"]<(-pi/2)){
		angle["v"]=-pi-angle["v"]
	}
	if(snout.pos["x"]<0){
		angle["h"]=pi+angle["h"]
	}
	eye.pos=head.pos+snout.pos*eye.frac
	return(list(angle=angle, eye.pos=eye.pos))
}
	
v.angle<-function(eye.loc, object, eye.angle){
	##determines the 3D-angle (in radians) between a viewing direction from an eye with a given location
		##to an object with a given location
	##eye.loc and object indicate the position of the eye and the object in a cartesian coordinate system
		##must be named vectors; names must be x, y, and z
		##x and y denote horizontal (E-W) and vertical (N-S) position, and z the elevation
	##eye.angle indicates the viewing direction of the dog (in radians)
		##must be named h and and v for the horizontal (x-y) and vertical (z) angle
		##h=0 means looking to the east; v<0 and >0 mean looking down- and upwards respectively
	if(!all(eye.loc==object)){
		object=object[c("x", "y", "z")]
		eye.loc=eye.loc[c("x", "y", "z")]
		object=object-eye.loc
		eye.loc[]=0
		dist.obj.eye=sqrt(sum(object^2))
		object=object/dist.obj.eye
		object.angle=c(h=as.vector(atan(abs(object["y"])/abs(object["x"]))), v=as.vector(atan(object["z"]/sqrt(sum(object[c("x", "y")]^2)))))
		#browser()
		if(is.nan(object.angle["h"])){
			object.angle["h"]=0
		}else if(object["x"]<0 & object["y"]>=0){
			object.angle["h"]=pi-object.angle["h"]
		}else if(object["x"]<0 & object["y"]<0){
			object.angle["h"]=pi+object.angle["h"]
		}else if(object["x"]>=0 & object["y"]<0){
			object.angle["h"]=-object.angle["h"]
		}
		#XX = sin(N1 * pi / 180) * sin(N2 * pi / 180) + cos(N1 * pi / 180) * cos(N2 * pi / 180) * cos(abs(E1 - E2) * pi / 180)
		##1->eye
		XX = sin(eye.angle["v"]) * sin(object.angle["v"]) + 
			cos(eye.angle["v"]) * cos(object.angle["v"]) * cos(abs(eye.angle["h"] - object.angle["h"]))
		if(eye.angle["v"] != object.angle["v"] | eye.angle["h"] != object.angle["h"]){
			XX = atan(-XX / sqrt(-XX * XX + 1)) + pi/2
		}else{
			XX = 0
		}
		return(XX)
	}else{
		return(NA)
	}
}


