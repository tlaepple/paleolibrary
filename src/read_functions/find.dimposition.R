

find.dimposition<-function(data.nc,varname,dimname)
  # find.dimposition
  # searches position of dimension in variable
{
  
   for (i in 1:length(data.nc$var))
   {
     if (data.nc$var[[i]]$name==varname)
        for (j in 1:length(data.nc$var[[i]]$dim))
       {
           if(is.element(data.nc$var[[i]]$dim[[j]]$name,dimname)) return(j)
       }
   }
   # if not found ...  
   return(NULL)
  
}
  
# debugging
# find.dimposition(ex.nc,"wisoaprt_d","time")