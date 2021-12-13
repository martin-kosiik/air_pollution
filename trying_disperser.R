devtools::install_github( 'lhenneman/SplitR')
devtools::install_github( 'lhenneman/disperseR')
devtools::install_github( 'lhenneman/hyspdisp')

library(hyspdisp)
library(SplitR)
library(disperseR)
'disperser_dir'
library(here)
setwd(here::here())

browseVignettes(package = "hyspdisp")


#install.packages("raster", version='3.4-5')

library(raster)
library(raster, version='3.4-5')
library(raster, lib.loc = 'C:/Users/marti/OneDrive/Dokumenty/R/win-library/4.1-alt')
library(rasteralt, lib.loc = 'C:/Users/marti/OneDrive/Dokumenty/R/win-library/4.1-alt')

detach("package:raster", unload=TRUE)

unloadNamespace("raster_3.5-2")

sessionInfo()

library(ncdf4)
my.file <- nc_open("C:\\Users\\marti\\Dropbox\\research_projects\\air_pollution\\disperser_dir\\main\\input\\hpbl\\hpbl.mon.mean.nc")


detachAllPackages <- function() {
  
  basic.packages <- c("package:stats","package:graphics","package:grDevices","package:utils","package:datasets","package:methods","package:base")
  
  package.list <- search()[ifelse(unlist(gregexpr("package:",search()))==1,TRUE,FALSE)]
  
  package.list <- setdiff(package.list,basic.packages)
  
  if (length(package.list)>0)  for (package in package.list) detach(package, character.only=TRUE)
  
}

detachAllPackages()


install.packages('raster')
disperseR::create_dirs(location = 'disperser_dir') 

require(remotes)
install_version("raster", version = "3.4-5", repos = "http://cran.us.r-project.org", lib = "C:/Users/marti/OneDrive/Dokumenty/R/win-library/4.1-alt")


library(raster)
r <- raster(ncol=10, nrow=10, xmx=-80, xmn=-150, ymn=20, ymx=60)
values(r) <- runif(ncell(r))

plot(r)

disperseR::get_data(data = "all", 
                    start.year = "2005", 
                    start.month = "11", 
                    end.year = "2006", 
                    end.month = "02")
View(disperseR::units)

pblheight <- disperseR::get_data(data = "pblheight")

hpbl_rasterin <- raster::brick(x = "C:\\Users\\marti\\Dropbox\\research_projects\\air_pollution\\disperser_dir\\main\\input\\hpbl\\air.mon.mean.nc")

plot(hpbl_rasterin)

f <- system.file("C:\\Users\\marti\\Dropbox\\research_projects\\air_pollution\\disperser_dir\\main\\input\\hpbl\\hpbl.mon.mean.nc", package="raster")
pack <- available.packages()
pack["disperseR","Depends"]

pack %>% as_tibble() %>% 
  filter(str_detect(Package, "disperseR"))

directory <- hpbl_dir
file_path <- file.path(directory, "hpbl.mon.mean.nc")
url <- "https://www.esrl.noaa.gov/psd/repository/entry/get/hpbl.mon.mean.nc?entryid=synth%3Ae570c8f9-ec09-4e89-93b4-babd5651e7a9%3AL05BUlIvTW9udGhsaWVzL21vbm9sZXZlbC9ocGJsLm1vbi5tZWFuLm5j"

Sys.setenv(TZ = "UTC")
hpbl_rasterin <- raster::brick(x = file_path, varname = "hpbl")
raster::crs(hpbl_rasterin) <- "+proj=lcc +x_0=5632642.22547 +y_0=4612545.65137 +lat_0=50 +lon_0=-107 +lat_1=50"



unitsrun2005 <- disperseR::units %>% 
  dplyr::filter(year == 2005) %>% # only get data for 2005
  dplyr::top_n(2, SOx)  # sort and take the two rows with the biggest value for SOx

unitsrun2006 <- disperseR::units %>% 
  dplyr::filter(year == 2006) %>%  # only get data for 2006
  dplyr::top_n(2, SOx)  # sort and take the two rows with the biggest value for SOx

head(unitsrun2005)
#>        ID Latitude Longitude      SOx      CO2      NOx   Height inputed year    uID
#> 1: 3136-1  40.6604  -79.3411 183555.8 13322458 13787.34 171.6024       0 2005 3136.1
#> 2: 3149-1  41.0714  -76.6672 240250.9 18188647 22388.27 213.3600       0 2005 3149.1

# append together and transform to data table 
unitsrun<-data.table::data.table(rbind(unitsrun2005, unitsrun2006))


input_refs <- disperseR::define_inputs(units = unitsrun,
                                       startday = '2005-11-01',
                                       endday = '2006-02-28',
                                       start.hours =  c(0, 6, 12, 18),
                                       duration = 120)

head(input_refs)


input_refs_subset <- input_refs[format(as.Date(input_refs$start_day,
                                               format = "%Y-%m-%d"),
                                       format = "%d") == "01" & start_hour == 0]


hysp_raw <- disperseR::run_disperser_parallel(input.refs = input_refs_subset,
                                              pbl.height = 'jnjn',
                                              overwrite = T,
                                              species = 'so2',
                                              proc_dir = proc_dir,
                                              npart = 100,
                                              keep.hysplit.files = FALSE, ## FALSE BY DEFAULT
                                              #mc.cores = parallel::detectCores()
                                              mc.cores = 1)


# trims particles that are above the global max boundary value
#disp_df_trim <- disp_df[height <= 2665]

