# Load Libraries
library(rjson)
library(tidycensus)
library(tidyverse)
library(sf)
library(spdep)
library(caret)
library(ckanr)
library(FNN)
library(grid)
library(gridExtra)
library(ggcorrplot)
library(jtools)  
library(viridis)
library(kableExtra)
library(rlist)
library(dplyr)
library(osmdata)
library(geosphere)
library(fastDummies)
options(scipen=999)
options(tigris_class = "sf")


# Themes and Functions
mapTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle=element_text(face="italic"),
    plot.caption=element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),axis.title = element_blank(),
    axis.text = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2)
  )
}

plotTheme <- function(base_size = 12) {
  theme(
    text = element_text( color = "black"),
    plot.title = element_text(size = 14,colour = "black"),
    plot.subtitle = element_text(face="italic"),
    plot.caption = element_text(hjust=0),
    axis.ticks = element_blank(),
    panel.background = element_blank(),
    panel.grid.major = element_line("grey80", size = 0.1),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA, size=2),
    strip.background = element_rect(fill = "grey80", color = "white"),
    strip.text = element_text(size=12),
    axis.title = element_text(size=12),
    axis.text = element_text(size=10),
    plot.background = element_blank(),
    legend.background = element_blank(),
    legend.title = element_text(colour = "black", face = "italic"),
    legend.text = element_text(colour = "black", face = "italic"),
    strip.text.x = element_text(size = 14)
  )
}

palette5 <- c("#25CB10", "#5AB60C", "#8FA108",   "#C48C04", "#FA7800")

qBr <- function(df, variable, rnd) {
  if (missing(rnd)) {
    as.character(quantile(round(df[[variable]],0),
                          c(.01,.2,.4,.6,.8), na.rm=T))
  } else if (rnd == FALSE | rnd == F) {
    as.character(formatC(quantile(df[[variable]]), digits = 3),
                 c(.01,.2,.4,.6,.8), na.rm=T)
  }
}

q5 <- function(variable) {as.factor(ntile(variable, 5))}

# Functions

nn_function <- function(measureFrom,measureTo,k) {
  measureFrom_Matrix <- as.matrix(measureFrom)
  measureTo_Matrix <- as.matrix(measureTo)
  nn <-   
    get.knnx(measureTo, measureFrom, k)$nn.dist
  output <-
    as.data.frame(nn) %>%
    rownames_to_column(var = "thisPoint") %>%
    gather(points, point_distance, V1:ncol(.)) %>%
    arrange(as.numeric(thisPoint)) %>%
    group_by(thisPoint) %>%
    summarize(pointDistance = mean(point_distance)) %>%
    arrange(as.numeric(thisPoint)) %>% 
    dplyr::select(-thisPoint) %>%
    pull()
  
  return(output)  
}

#Function MultipleRingBuffer
multipleRingBuffer <- function(inputPolygon, maxDistance, interval) 
{
  #create a list of distances that we'll iterate through to create each ring
  distances <- seq(0, maxDistance, interval)
  #we'll start with the second value in that list - the first is '0'
  distancesCounter <- 2
  #total number of rings we're going to create
  numberOfRings <- floor(maxDistance / interval)
  #a counter of number of rings
  numberOfRingsCounter <- 1
  #initialize an otuput data frame (that is not an sf)
  allRings <- data.frame()
  
  #while number of rings  counteris less than the specified nubmer of rings
  while (numberOfRingsCounter <= numberOfRings) 
  {
    #if we're interested in a negative buffer and this is the first buffer
    #(ie. not distance = '0' in the distances list)
    if(distances[distancesCounter] < 0 & distancesCounter == 2)
    {
      #buffer the input by the first distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #different that buffer from the input polygon to get the first ring
      buffer1_ <- st_difference(inputPolygon, buffer1)
      #cast this sf as a polygon geometry type
      thisRing <- st_cast(buffer1_, "POLYGON")
      #take the last column which is 'geometry'
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add a new field, 'distance' so we know how far the distance is for a give ring
      thisRing$distance <- distances[distancesCounter]
    }
    
    
    #otherwise, if this is the second or more ring (and a negative buffer)
    else if(distances[distancesCounter] < 0 & distancesCounter > 2) 
    {
      #buffer by a specific distance
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create the next smallest buffer
      buffer2 <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #This can then be used to difference out a buffer running from 660 to 1320
      #This works because differencing 1320ft by 660ft = a buffer between 660 & 1320.
      #bc the area after 660ft in buffer2 = NA.
      thisRing <- st_difference(buffer2,buffer1)
      #cast as apolygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #get the last field
      thisRing <- as.data.frame(thisRing$geometry)
      #create the distance field
      thisRing$distance <- distances[distancesCounter]
    }
    
    #Otherwise, if its a positive buffer
    else 
    {
      #Create a positive buffer
      buffer1 <- st_buffer(inputPolygon, distances[distancesCounter])
      #create a positive buffer that is one distance smaller. So if its the first buffer
      #distance, buffer1_ will = 0. 
      buffer1_ <- st_buffer(inputPolygon, distances[distancesCounter-1])
      #difference the two buffers
      thisRing <- st_difference(buffer1,buffer1_)
      #cast as a polygon
      thisRing <- st_cast(thisRing, "POLYGON")
      #geometry column as a data frame
      thisRing <- as.data.frame(thisRing[,ncol(thisRing)])
      #add teh distance
      thisRing$distance <- distances[distancesCounter]
    }  
    
    #rbind this ring to the rest of the rings
    allRings <- rbind(allRings, thisRing)
    #iterate the distance counter
    distancesCounter <- distancesCounter + 1
    #iterate the number of rings counter
    numberOfRingsCounter <- numberOfRingsCounter + 1
  }
  
  #convert the allRings data frame to an sf data frame
  allRings <- st_as_sf(allRings)
}

## Read Miami Data

MiamiDF <- st_read("studentsData.geojson")

MiamiDF <- st_transform(MiamiDF,'ESRI:102658')

MiamiDF <- dplyr::select(MiamiDF,-saleQual,-WVDB,-HEX,-GPAR,-County.2nd.HEX, 
                         -County.Senior,-County.LongTermSenior,-County.Other.Exempt,
                         -City.2nd.HEX,-City.Senior,-City.LongTermSenior,
                         -City.Other.Exempt,-MillCode,-Land.Use,
                         -Owner1,-Owner2,-Mailing.Address,-Mailing.City,
                         -Mailing.State,-Mailing.Zip,-Mailing.Country,
                         -Year,-Land,-Bldg,-Total,-Assessed,-County.Taxable,
                         -City.Taxable,-Legal1,-Legal2,-Legal3,-Legal4,-Legal5,-Legal6,-Units)

# Wrangle XF Columns
MiamiDF<- mutate(MiamiDF,XFs=paste(XF1,XF2,XF3, sep=",")) %>%
  dummy_cols(select_columns="XFs", split=",")

MiamiDF$LuxuryPool<- ifelse(MiamiDF$`XFs_Luxury Pool - Best`==1 | MiamiDF$`XFs_Luxury Pool - Better`==1 | MiamiDF$`XFs_Luxury Pool - Good.`==1,1,0)

MiamiDF$'3to6ftPool' <- ifelse(MiamiDF$`XFs_Pool COMM AVG 3-6' dpth`==1|MiamiDF$`XFs_Pool COMM BETTER 3-6' dpth`==1,1,0)

MiamiDF$'3to8ftPool' <- ifelse(MiamiDF$`XFs_Pool 6' res BETTER 3-8' dpth`==1|MiamiDF$`XFs_Pool 6' res AVG 3-8' dpth`==1,1,0)

MiamiDF <- rename(MiamiDF,`8ftres3to8ftPool`=`XFs_Pool 8' res BETTER 3-8' dpth`)

MiamiDF <- rename(MiamiDF, Whirpool=`XFs_Whirlpool - Attached to Pool (whirlpool area only)`)

MiamiDF <- rename(MiamiDF,`2to4ftPool`= `XFs_Pool - Wading - 2-4' depth`)

MiamiDF <- dplyr::select(MiamiDF,-"XFs_Pool 6' res BETTER 3-8' dpth",
                         -"XFs_Pool 6' res AVG 3-8' dpth",-"XFs_Luxury Pool - Best",
                         -"XFs_Luxury Pool - Better",-"XFs_Luxury Pool - Good.",
                         -"XFs_Pool COMM AVG 3-6' dpth",-"XFs_Pool COMM BETTER 3-6' dpth",
                         -"XFs_large",-"XFs_elec",-"XFs_plumb",-"XFs_Tiki Hut - Standard Thatch roof w/poles",
                         -"XFs_Tiki Hut - Good Thatch roof w/poles & electric",-"XFs_Tiki Hut - Better Thatch roof",
                         -"XFs_Bomb Shelter - Concrete Block",-"XFs_Tennis Court - Asphalt or Clay" ) 


# Wrangle Fence Columns
MiamiDF$Fence <- ifelse(MiamiDF$`XFs_Aluminum Modular Fence`==1| MiamiDF$`XFs_Wood Fence`==1|MiamiDF$`XFs_Chain-link Fence 4-5 ft high`==1|MiamiDF$`XFs_Wrought Iron Fence`==1|MiamiDF$`XFs_Chain-link Fence 6-7 ft high`==1|MiamiDF$`XFs_Concrete Louver Fence`==1|MiamiDF$`XFs_Chain-link Fence 8-9 ft high`==1|MiamiDF$`XFs_Glass fences in backyard applications`==1,1,0)

MiamiDF <- dplyr::select(MiamiDF,-"XFs_Aluminum Modular Fence",
                         -"XFs_Wood Fence",-"XFs_Chain-link Fence 4-5 ft high",
                         -"XFs_Wrought Iron Fence",-"XFs_Chain-link Fence 6-7 ft high",
                         -"XFs_Concrete Louver Fence",-"XFs_Chain-link Fence 8-9 ft high",
                         -"XFs_Glass fences in backyard applications")

# Wrangle Patio Columns
MiamiDF$Patio <- ifelse(MiamiDF$`XFs_Patio - Concrete Slab`==1|MiamiDF$`XFs_Patio - Concrete Slab w/Roof Aluminum or Fiber`==1|MiamiDF$`XFs_Patio - Wood Deck`==1|MiamiDF$`XFs_Patio - Marble`==1|MiamiDF$`XFs_Patio - Terrazzo`==1|MiamiDF$`XFs_Patio - Screened over Concrete Slab`==1|MiamiDF$`XFs_Patio - Exotic hardwood`==1|MiamiDF$`XFs_Patio - Concrete stamped or stained`==1,1,0)

MiamiDF <- dplyr::select(MiamiDF,-"XFs_Patio - Concrete Slab",-"XFs_Patio - Concrete Slab w/Roof Aluminum or Fiber",
                         -"XFs_Patio - Wood Deck",-"XFs_Patio - Marble",
                         -"XFs_Patio - Terrazzo",-"XFs_Patio - Screened over Concrete Slab",
                         -"XFs_Patio - Exotic hardwood",-"XFs_Patio - Concrete stamped or stained")

## Wrangle Dock Columns
MiamiDF$Docks <- ifelse(MiamiDF$`XFs_Dock - Wood on Light Posts`==1|MiamiDF$`XFs_Loading Dock/ Platform`==1|MiamiDF$`XFs_Dock - Concrete Griders on Concrete Pilings`==1|MiamiDF$`XFs_Dock - Wood Girders on Concrete Pilings`==1|MiamiDF$`XFs_Dock - Steel Pilings`==1,1,0)