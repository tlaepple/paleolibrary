#Function to interpolate a field to a given point
#It uses the nearest neighbour if the adjancents points are missing, if not it uses bilinear interpolation
#Experimental phase, only for 2D fields in the moment
#neue Funktion, erst in 2D konvertieren, dann verlängern um NA Probleme zu vermeiden

library(fUtilities)

selspace.interpolate<-function(data,lat1,lon1,SBOX=5)
{
  

  	  temp<-attributes(data)

	  if (dim(data)[1]>1) stop("Works only for 2D Fields (= only one timestep)")
        if (prod(dim(data)) != length(temp$lon)*length(temp$lat)) stop("N(data) != N(lat)*N(lon)")

       data2d<-matrix(data,length(temp$lon),length(temp$lat)) #make a 2D array
        
        #arrange to get a continous field
        d<-diff(temp$lon)
        if (max(d) > (min(d)+0.01))
          { nlon<-length(temp$lon)
            edgelon<-which(d==max(d))
            data2d<-rbind(data2d[(edgelon+1):nlon,],data2d[1:edgelon,])
            temp$lon<-c(temp$lon[(edgelon+1):nlon],temp$lon[1:edgelon]+360)
          }
        
        ### Copy the data 3 times....

          nlon<-length(temp$lon)
          data2d<-rbind(data2d[1:nlon,],data2d[1:nlon,],data2d[1:nlon,])
          temp$lon<-c(temp$lon[1:nlon]-360,temp$lon[1:nlon],temp$lon[1:nlon]+360)
        
        if (temp$lat[2] < temp$lat[1])    #if the latitudes are from + to -, reverse them
         {
          temp$lat<-rev(temp$lat)
          data2d<-data2d[,rev(seq(len=ncol(data2d)))]
        }
   



    #attention... midpoints are given...


    indexLat<-rev(which(temp$lat<=lat1))[1]:which(temp$lat>=lat1)[1]
    indexLon<-rev(which(temp$lon<=lon1))[1]:which(temp$lon>=lon1)[1]
 
    #First check if we have any missing data
    data.select<-data2d[indexLon,indexLat]

    if (sum(is.na(data.select)) > 0)
	{ 
					
			#Create an area around the point
			tempindexLat<-((indexLat[1]-SBOX):(indexLat[1]+SBOX))
		      tempindexLon<-((indexLon[1]-SBOX):(indexLon[1]+SBOX))

			#remove areas outside the boundaries
			tempindexLat<-tempindexLat[tempindexLat>0]
			tempindexLon<-tempindexLon[tempindexLon>0]
		      tempindexLat<-tempindexLat[tempindexLat<=length(temp$lat)]
			tempindexLon<-tempindexLon[tempindexLon<=length(temp$lon)]

			#Get the nearest neighbours
			res<-expand.grid(tempindexLon,tempindexLat)

			#only retain the nearest nonmissing neighbours
			indexV<-!is.na(diag(data2d[res[,1],res[,2]]))
			
			if (sum(indexV) == 0) return(NA)
			x<-res[indexV,1]
			y<-res[indexV,2]

			D2i <- (temp$lon[x] - lon1)^2 + (temp$lat[y] - lat1)^2  #Distances to the points in the area
                                 
			neighbour <- order(D2i)[1]  #nearest neighbour is the one with the smallest distance
			result<-data2d[x[neighbour],y[neighbour]]

	}

		


	 else
{

    if ((length(indexLat) == 1) & (length(indexLon) == 1)) result<-data.select else #we are on the point, no interpolation
    if (length(indexLat) == 1)  result<-(approx(temp$lon[indexLon],data.select,lon1)$y) else #we are on a horizontal line 1D Problem
    if (length(indexLon) == 1)  result<-(approx(temp$lat[indexLat],data.select,lat1)$y) else #we are on a vertical line, 1D Problem
   {
	g<-  expand.grid(temp$lon[indexLon],temp$lat[indexLat])  #Normal case, 2D problem
	result<-linearInterpp(g[,1],g[,2],data.select,lon1,lat1)$z
   }
}


return(result)
}

