# Land change modeling based on Weight of Evidence for spatial distribution
# steps: 
        #1 transition matrix based on LULC, 
        #2 WoE based on variables, 
        #3 predict (prob and LU allocation)

# load packages and data ####
  library(raster)
  library(sp)
  library(rgdal)
  library(gstat)
    
  library(dplyr)
  
  lulc00 = raster("data/sample_mapb2000.tif")
    rcl <- matrix(c(3,1, 13,1, 15,2, 19,2, 24,3, 33, 4), ncol=2, byrow=TRUE) #1=veg, 2=agric, 3=urban, 4=water
    lulc00 <- reclassify(lulc00, rcl) #reclassify the LULC to simply the operations
  lulc18 = raster("data/sample_mapb2018.tif")
    lulc18 <- reclassify(lulc18, rcl)
    remove(rcl)
    
  DEM = raster("data/dem.tif")
  lots = readOGR("data/lots.shp")
  roads = readOGR("data/main_road.shp")
  soils = readOGR("data/solis_utm.shp")
    rasterize(soils, DEM, field=as.factor(soils$SOIL_ID), filename="data/terrain", overwrite=T)
    terrain = raster("data/terrain.grd")

  # distance raster
      rcl <- matrix(c(1,NA, 2,2, 3,3, 4,4), ncol=2, byrow=TRUE)
      nfRaster <- reclassify(lulc00, rcl) #reclassify the LULC to simplify the operations
        remove(rcl)
    agrDist = distance(nfRaster)
    
      rcl <- matrix(c(1,NA, 2,NA, 3,3, 4,NA), ncol=2, byrow=TRUE)
      uRaster <- reclassify(lulc00, rcl) #reclassify the LULC to simplify the operations
        remove(rcl)
    urbanDist = distance(uRaster)
    
    roadDist = rasterize(roads, DEM)
    roadDist = distance(roadDist)
    
    
  # Rasters as data frame
    s = stack(list(lulc00=lulc00, lulc18=lulc18, DEM=DEM, terrain=terrain,
                   roadDist=roadDist, agrDist=agrDist, urbanDist=urbanDist))
    s = as.data.frame(s, xy=TRUE)
      s=na.omit(s)
    #rasterToPoints(s) # back to raster
 
# 1. Transition matrix/probabilities ####
  tm = matrix(table(s[,4:3]),ncol = 4) # transition matrix; note in 's[,4:3]' that R read from column to row
  tp = tm/rowSums(tm)
  
  # 1.1. Markov chain ####
    # see the visualizatoin script to depict markov chain probabilities
    # Also, the df 'steps' and the 'initiate probability vectors' has the values to be used for the projected risk using WoE

# 2. WoE and land change (LC) probabilities ####
  library(Information)
  #library(gridExtra)
    
  # create a binary output, in order to calculate the WoE
    s_def = subset(s, lulc00==1)
      s_def$def = s_def$lulc18-s_def$lulc00
      s_def$def[s_def$def>1] = 1
  
  # WoE and posterior prob
  IV = create_infotables(data=s_def, y="def", bins=5) # needs to incorporate binary variables as section #4
    IV$Summary #see output; IV is relevant to ranking variables' weight
      #print(IV$Tables$DEM, row.names=FALSE)
      
    # 2.1. Posterior probabilities of LC, from WoE ####
      # from veg to agric according to a given variable (i.e., DEM)
      # save into WOE results to keep the bins of each variable
      IV$Tables$DEM$pDEM = exp(IV$Tables$DEM[,4]+log(tp[2,1]/(1-tp[2,1])))/
                                (1+exp(IV$Tables$DEM[,4]+log(tp[2,1]/(1-tp[2,1]))))
      IV$Tables$terrain$pterrain = exp(IV$Tables$terrain[,4]+log(tp[2,1]/(1-tp[2,1])))/
                                        (1+exp(IV$Tables$terrain[,4]+log(tp[2,1]/(1-tp[2,1]))))
      IV$Tables$roadDist$proadDist = exp(IV$Tables$roadDist[,4]+log(tp[2,1]/(1-tp[2,1])))/
                                          (1+exp(IV$Tables$roadDist[,4]+log(tp[2,1]/(1-tp[2,1]))))
      IV$Tables$agrDist$pagrDist = exp(IV$Tables$agrDist[,4]+log(tp[2,1]/(1-tp[2,1])))/
                                          (1+exp(IV$Tables$agrDist[,4]+log(tp[2,1]/(1-tp[2,1]))))
      IV$Tables$urbanDist$puDist = exp(IV$Tables$urbanDist[,4]+log(tp[2,1]/(1-tp[2,1])))/
                                        (1+exp(IV$Tables$urbanDist[,4]+log(tp[2,1]/(1-tp[2,1]))))
      
    # 2.2. Transfer the probabilities to the stack of rasters ####
      library(readr)
      s$pDEM = ifelse(s$DEM<parse_number(IV$Tables$DEM[2,1], locale = locale(grouping_mark = ";")), IV$Tables$DEM[1,6],
                      ifelse(s$DEM<parse_number(IV$Tables$DEM[3,1], locale = locale(grouping_mark = ";")), IV$Tables$DEM[2,6],
                      ifelse(s$DEM<parse_number(IV$Tables$DEM[4,1], locale = locale(grouping_mark = ";")), IV$Tables$DEM[3,6],
                      ifelse(s$DEM<parse_number(IV$Tables$DEM[5,1], locale = locale(grouping_mark = ";")), IV$Tables$DEM[4,6],
                      ifelse(s$DEM<352, IV$Tables$DEM[5,6],NA
                        )))))
      s$pterrain = ifelse(s$terrain==1, IV$Tables$terrain[1,6],
                          ifelse(s$terrain<4, IV$Tables$terrain[2,6],
                          ifelse(s$terrain<5, IV$Tables$terrain[3,6],
                          ifelse(s$terrain<7, IV$Tables$terrain[4,6], NA
                            ))))
      s$pRoad = ifelse(s$roadDist<parse_number(IV$Tables$roadDist[2,1], locale = locale(grouping_mark = ";")), IV$Tables$roadDist[1,6],
                          ifelse(s$roadDist<parse_number(IV$Tables$roadDist[3,1], locale = locale(grouping_mark = ";")), IV$Tables$roadDist[2,6],
                          ifelse(s$roadDist<parse_number(IV$Tables$roadDist[4,1], locale = locale(grouping_mark = ";")), IV$Tables$roadDist[3,6],
                          ifelse(s$roadDist<parse_number(IV$Tables$roadDist[5,1], locale = locale(grouping_mark = ";")), IV$Tables$roadDist[4,6], 
                          ifelse(s$roadDist<13000, IV$Tables$roadDist[5,6], NA
                                        )))))
      s$pAgr = ifelse(s$agrDist<parse_number(IV$Tables$agrDist[2,1], locale = locale(grouping_mark = ";")), IV$Tables$agrDist[1,6],
                         ifelse(s$agrDist<parse_number(IV$Tables$agrDist[3,1], locale = locale(grouping_mark = ";")), IV$Tables$agrDist[2,6],
                         ifelse(s$agrDist<parse_number(IV$Tables$agrDist[4,1], locale = locale(grouping_mark = ";")), IV$Tables$agrDist[3,6],
                         ifelse(s$agrDist<parse_number(IV$Tables$agrDist[5,1], locale = locale(grouping_mark = ";")), IV$Tables$agrDist[4,6],
                         ifelse(s$agrDist<1400, IV$Tables$agrDist[5,6],NA
                                )))))
      s$pUrb = ifelse(s$urbanDist<parse_number(IV$Tables$urbanDist[2,1], locale = locale(grouping_mark = ";")), IV$Tables$urbanDist[1,6],
                      ifelse(s$urbanDist<parse_number(IV$Tables$urbanDist[3,1], locale = locale(grouping_mark = ";")), IV$Tables$urbanDist[2,6],
                      ifelse(s$urbanDist<parse_number(IV$Tables$urbanDist[4,1], locale = locale(grouping_mark = ";")), IV$Tables$urbanDist[3,6],
                      ifelse(s$urbanDist<parse_number(IV$Tables$urbanDist[5,1], locale = locale(grouping_mark = ";")), IV$Tables$urbanDist[4,6],
                      ifelse(s$urbanDist<12000, IV$Tables$agrDist[5,6],NA
                                           )))))
      
      s$prob = ifelse(s$lulc00==1, rowMeans(s[,10:length(s)], na.rm=TRUE), NA)
        s = na.omit(s)
        
      # Crete a map (raster) of posterior probabilities
      Dp = rasterize(s[,1:2], terrain, field=s$prob, filename="data/Dprob", overwrite=T)
      
# 3. LULC allocation (predict) ####
  # for the LU allocation, combine WoE and transition matrix probabilities
    # i.e., gross rate of def is rowSums(tm)[1]*tp[2,1] = 39836 pts from 2000 to 2018
      luc = s %>%
        slice_max(prob,n=rowSums(tm)[1]*tp[2,1], with_ties=F)  %>%
        mutate(luc=1)
      
      lchange = rasterize(luc[,1:2], terrain, field=luc$luc)
        plot(lchange, col="Red")
      
      # checking real def occurrence
        c_def = subset(s, lulc00==1 & lulc18==2)
          c_def$def = c_def$lulc18-c_def$lulc00 #just for checking with: table(c_def$def)
        def_pos = rasterize(c_def[,1:2], terrain, field=c_def$def)
          plot(def_pos, col="Red")
      
      # matching def rate and spatial distribution
      count(luc)/count(c_def) # comparet our estimate to actual deforestation
      crosstab(lchange, def_pos)/count(c_def) # accuracy in the distribution, but doesn't count false-positives
      
      # Replicate for regrowth
      
      
# 4. Tests with logit ####
    # Call the binary df from WoE section; it helps with thr R2 for this estimates (i.e., 10%, which is not the best)
      reg = lm(def ~DEM +as.factor(terrain) +agrDist +roadDist, data = s_def, family = binomial)
      summary(reg)
  