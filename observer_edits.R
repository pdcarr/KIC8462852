######## removing obvious wild points user by user
editUser <- data.frame(obsCode= "UJHA",startJD=2457650,endJD=2457760,band="V",stringsAsFactors=FALSE)
editUser <- rbind(editUser,c("PXR",2457512,2457513,"V"),c("JSJA",2457646,2457647,"V"),c("SGEA",2457879,2457880,"V"),c("SDB",2457609,2457610,"B"))
editUser <- rbind(editUser,c("DUBF",2457737,2457739,"B"),c("DUBF",2457737,2457739,"V"),c("DUBF",2457737.22213,2457739.22213,"B"))
