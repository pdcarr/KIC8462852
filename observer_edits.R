######## removing obvious wild points user by user
editUser <- data.frame(obsCode= "UJHA",startJD=2457650,endJD=2457760,band="V",stringsAsFactors=FALSE)
editUser <- rbind(editUser,c("PXR",2457512,2457513,"V"),c("JSJA",2457646,2457647,"V"),c("SGEA",2457879,2457880,"V"),c("SDB",2457609,2457610,"B"))
editUser <- rbind(editUser,c("DUBF",2457737,2457739,"B"),c("DUBF",2457737,2457739,"V"),c("DUBF",2457737.22213,2457739.22213,"B"))
editUser <- rbind(editUser,c("SJAR", 2457620.8,2457621.5,"R"),c("DUBF", 2457737.5, 2457738.8,"R"),c("LDJ",2457392.44315,2457394.44315,"V"))
editUser <- rbind(editUser,c("JM",2457466.9,2457468.6,"V"),c("PXR",2457573.0,2457574.5,"V"),c("MJB",2457634.0,2457635.3,"V"))
editUser <- rbind(editUser,c("HDHA",2457907,2457908.1,"V"),c("UJHA",2457925.2,2457926.5,"V"))