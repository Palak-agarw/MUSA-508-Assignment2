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

# Wrangle Miami Data
MiamiDF <- st_read("studentsData.geojson")

MiamiDF <- st_transform(MiamiDF,'ESRI:102658')

MiamiDF <- dplyr::select(MiamiDF,-saleQual,-WVDB,-HEX,-GPAR,-County.2nd.HEX,
                         -County.Senior,-County.LongTermSenior,-County.Other.Exempt,
                         -City.2nd.HEX,-City.Senior,-City.LongTermSenior,
                         -City.Other.Exempt,-MillCode,-Land.Use,
                         -Owner1,-Owner2,-Mailing.Address,-Mailing.City,
                         -Mailing.State,-Mailing.Zip,-Mailing.Country)

# Wrangle XF Columns
MiamiDF<- mutate(MiamiDF,XFs=paste(XF1,XF2,XF3, sep=",")) %>%
  dummy_cols(select_columns="XFs", split=",")

MiamiDF <- st_as_sf(MiamiDF)

# Join Neighborhood Data
Neighborhoods <- st_read("https://opendata.arcgis.com/datasets/2f54a0cbd67046f2bd100fb735176e6c_0.geojson")

Neighborhoods <- st_transform(Neighborhoods,'ESRI:102658')

Municipality <- st_read('https://opendata.arcgis.com/datasets/bd523e71861749959a7f12c9d0388d1c_0.geojson')

Municipality <- st_transform(Municipality,'ESRI:102658')

MiamiDF <- st_join(MiamiDF, Neighborhoods, join = st_intersects) 

# Remove Challenge Houses

MiamiDFKnown <- MiamiDF[!(MiamiDF$SalePrice==0),]

# Map House Sales and Neighborhoods

ggplot() +
  geom_sf(data=Municipality)+
  geom_sf(data=MiamiDF)


ggplot() +
  geom_sf(data = Neighborhoods, fill = "grey40") +
  geom_sf(data = MiamiDFKnown, aes(colour = q5(SalePrice)), 
          show.legend = "point", size = 1) +
  scale_colour_manual(values=palette5,
                      labels=qBr(MiamiDFKnown,"SalePrice"),
                      name="Quintile\nBreaks") +
  labs(title="Sale Price, Miami") +
  mapTheme()

#Transit Data 

metroStops <- st_read("https://opendata.arcgis.com/datasets/ee3e2c45427e4c85b751d8ad57dd7b16_0.geojson") 
metroStops <- metroStops %>% st_transform('ESRI:102658')

# Plot of the metro stops
ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data=MiamiDFKnown)+
  geom_sf(data=metroStops, 
          aes(colour = 'red' ),
          show.legend = "point", size= 1.2)

metroBuffers <- 
  rbind(
    st_buffer(metroStops, 2640) %>%
      mutate(Legend = "Buffer") %>%
      dplyr::select(Legend),
    st_union(st_buffer(metroStops, 2640)) %>%
      st_sf() %>%
      mutate(Legend = "Unioned Buffer"))

ggplot() +
  geom_sf(data=metroBuffers) +
  geom_sf(data=metroStops, show.legend = "point") +
  facet_wrap(~Legend) + 
  mapTheme()

# Create an sf object with ONLY the unioned buffer
buffer <- filter(metroBuffers, Legend=="Unioned Buffer")
buffer <- buffer %>% st_transform('ESRI:102658')

# Clip the Miami training DF ... by seeing which tracts intersect (st_intersection)
# with the buffer and clipping out only those areas
clip <- 
  st_intersection(buffer, MiamiDFKnown) %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Clip")

ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data = clip) 

# Do a spatial selection to see which tracts touch the buffer
selection <- 
  MiamiDFKnown[buffer,] %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Spatial Selection")

ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data = selection) +
  geom_sf(data = buffer, fill = "transparent", color = "red")

selectCentroids <-
  st_centroid(MiamiDFKnown)[buffer,] %>%
  st_drop_geometry() %>%
  left_join(dplyr::select(MiamiDFKnown, Folio)) %>%
  st_sf() %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Select by Centroids")

ggplot() +
  geom_sf(data = selectCentroids) +
  geom_sf(data = buffer, fill = "transparent", color = "red")+
  theme(plot.title = element_text(size=22))

MiamiDFKnown <- 
  rbind(
    st_centroid(MiamiDFKnown)[buffer,] %>%
      st_drop_geometry() %>%
      left_join(MiamiDFKnown) %>%
      st_sf() %>%
      mutate(TOD = 1),
    st_centroid(MiamiDFKnown)[buffer, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(MiamiDFKnown) %>%
      st_sf() %>%
      mutate(TOD = 0))

# Multiple Ring method

#Creating the buffer around the transit stops
MiamiDFKnown <-
  st_join(st_centroid(MiamiDFKnown), 
          multipleRingBuffer(st_union(metroStops), 47520, 1320)) %>%
  st_drop_geometry() %>%
  left_join(MiamiDFKnown) %>%
  st_sf() %>%
  mutate(Distance = distance / 5280)#convert to miles

MiamiDFKnown <- 
  MiamiDFKnown %>%
  mutate(NewDistance.cat = case_when(
    Distance >= 0 & Distance < 0.25  ~ "Quater Mile",
    Distance >= 0.25 & Distance < 0.5  ~ "Half Mile",
    Distance >= 0.5 & Distance < 0.75  ~ "Three Quater Mile",
    Distance > 1                    ~ "More than one Mile"))

# Correlation Matrix

numericVars <- 
  select_if(st_drop_geometry(MiamiDFKnown), is.numeric) %>% na.omit() %>%
  select(SalePrice.x,Bed, Bath, Stories, Units, YearBuilt, LivingSqFt, NewDistance.cat)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 

# Regression

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat))
summ(reg)
summary(reg)


# Bar/ restaurant data

miami.base <- 
  st_read("https://opendata.arcgis.com/datasets/5ece0745e24b4617a49f2e098df8117f_0.geojson") %>%
  filter(NAME == "MIAMI BEACH" | NAME == "MIAMI") %>%
  st_union()

xmin = st_bbox(miami.base)[[1]]
ymin = st_bbox(miami.base)[[2]]
xmax = st_bbox(miami.base)[[3]]  
ymax = st_bbox(miami.base)[[4]]

ggplot() +
  geom_sf(data=miami.base, fill="black") +
  geom_sf(data=st_as_sfc(st_bbox(miami.base)), colour="red", fill=NA) 

bars <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature(key = 'amenity', value = c("bar", "pub", "restaurant")) %>%
  osmdata_sf()

bars <- 
  bars$osm_points %>%
  .[miami.base,]

bars <- bars %>% st_transform('ESRI:102658')

ggplot() +
  geom_sf(data=miami.base, fill="black") +
  geom_sf(data=bars, colour="red", size=.75)

MiamiDFKnown <- MiamiDFKnown %>%
  mutate(geomtry = st_centroid(MiamiDFKnown))

## Nearest Neighbor Feature
st_c <- st_coordinates

MiamiDFKnown <-
  MiamiDFKnown %>% 
  mutate(
    bar_nn1 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(bars), 1),
    bar_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(bars), 5), 
    bar_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(bars), 10), 
    bar_nn4 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(bars), 15), 
    bar_nn5 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(bars), 20)) 

## bar cor
st_drop_geometry(MiamiTrainingDF) %>% 
  mutate(Age = 2020 - YearBuilt) %>%
  dplyr::select(SalePrice, starts_with("bar_")) %>%
  filter(SalePrice <= 1000000) %>%
  gather(Variable, Value, -SalePrice) %>% 
  ggplot(aes(Value, SalePrice)) +
  geom_point(size = .5) + geom_smooth(method = "lm", se=F, colour = "#FA7800") +
  facet_wrap(~Variable, nrow = 1, scales = "free") +
  labs(title = "Price as a function of continuous variables") +
  plotTheme()


# Correlation Matrix

numericVars <- 
  select_if(st_drop_geometry(MiamiTrainingDF), is.numeric) %>% na.omit() %>%
  select(SalePrice,Bed, Bath, Stories, Units, YearBuilt, LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 

# Regression
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5))
summ(reg)
summary(reg)

## Shoreline
Shoreline <- st_read("https://opendata.arcgis.com/datasets/58386199cc234518822e5f34f65eb713_0.geojson")
Shoreline<- st_transform(Shoreline,'ESRI:102658')

ggplot() +
  geom_sf(data = Neighborhoods, fill = "grey40") +
  geom_sf(data = Shoreline, color = "blue")

MiamiDFKnown <- MiamiDFKnown %>%
  mutate(distancetoshore= (st_distance(st_centroid(MiamiDFKnown), Shoreline) / 5280))

MiamiDFKnown <- MiamiDFKnown %>%
  mutate(distancetoshore= dist.mat <- geosphere::dist2Line(p = Shoreline , line = MiamiDFKnown))


MiamiDFKnown <- units::drop_units(MiamiDFKnown)

MiamiDFKnown <- MiamiDFKnown %>%
  mutate(Distancetoshore= as.numeric(distancetoshore))


#Coastline Data
Coastline<-opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature("natural", "coastline") %>%
  osmdata_sf()

MiamiDFKnown <- st_transform(MiamiDFKnown,'ESRI:37001' )

#add to MiamiDFKnown and convert to miles
MiamiDFKnown <-
  MiamiDFKnown %>%  
  mutate(CoastDist=(geosphere::dist2Line(p=st_coordinates(st_centroid(MiamiDFKnown)),
                                         line=st_coordinates(Coastline$osm_lines)[,1:2])*0.00062137)[,1])

MiamiDFKnown <- st_transform(MiamiDFKnown,'ESRI:102658' )

## Something weird going on with this regression
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
                     dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, 
                                   bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5, CoastDist))
summ(reg)
summary(reg)

# School Attendance Areas
## Elementary Schools
## Doesn't include Miami Beach
## Something is wrong with projection
ElementarySchool <- st_read("https://opendata.arcgis.com/datasets/19f5d8dcd9714e6fbd9043ac7a50c6f6_0.geojson")

ElementarySchool<- st_transform(ElementarySchool,'ESRI:102658') %>%
  select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-FLAG,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

ElementarySchool<-filter(ElementarySchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = ElementarySchool, fill = "grey40") +
  geom_sf(data = MiamiDF)+
  geom_sf(data = Neighborhoods)

## Middle Schools
MiddleSchool <- st_read("https://opendata.arcgis.com/datasets/dd2719ff6105463187197165a9c8dd5c_0.geojson")

MiddleSchool<- st_transform(MiddleSchool,'ESRI:102658') %>%
  select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

MiddleSchool<-filter(MiddleSchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = MiddleSchool, fill = "grey40") +
  geom_sf(data = MiamiDF)


## High Schools
HighSchool <- st_read("https://opendata.arcgis.com/datasets/9004dbf5f7f645d493bfb6b875a43dc1_0.geojson")

HighSchool<- st_transform(HighSchool,'ESRI:102658') %>%
  select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

HighSchool<-filter(HighSchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = HighSchool, fill = "grey40") +
  geom_sf(data = MiamiDF) +
  geom_sf(data = Neighborhoods)

ElementarySchool['Elementary'] <- "Elementary"
MiddleSchool['Middle'] <- "Middle"
HighSchool['High']<-"High"

ElementarySchool <- ElementarySchool %>% 
  select(-NAME,-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area,-geometry)

MiddleSchool <- MiddleSchool %>% 
  select(-NAME,-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area,-geometry)

HighSchool <- HighSchool %>% 
  select(-NAME,-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area,-geometry)

MiamiDFKnown <- st_join(MiamiDFKnown, ElementarySchool, join = st_intersects) 
MiamiDFKnown <- st_join(MiamiDFKnown, MiddleSchool, join = st_intersects) 
MiamiDFKnown <- st_join(MiamiDFKnown, HighSchool, join = st_intersects) 

MiamiDFKnown['AllSchool'] <-  paste(MiamiDFKnown$ELemantary, MiamiDFKnown$Middle, MiamiDFKnown$High, sep =',')

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
                     dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                                   bar_nn5, CoastDist, AllSchool))
summ(reg)
summary(reg)

## Parks 

Parks <- st_read("https://opendata.arcgis.com/datasets/8c9528d3e1824db3b14ed53188a46291_0.geojson")

Parks <- st_transform(Parks,'ESRI:102658')

MiamiDFKnown <-
  MiamiDFKnown %>% 
  mutate(
    park_nn1 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(Parks), 1),
    park_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(Parks), 3), 
    park_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(Parks), 4), 
    park_nn4 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(Parks), 5), 
    park_nn5 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(Parks), 10)) 

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, AllSchool, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5 ))
summ(reg)
summary(reg)

## Place of worship

worship <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature(key = 'amenity', value = c("place_of_worship")) %>%
  osmdata_sf()

worship <- 
  worship$osm_points %>%
  .[miami.base,]

worship <- worship %>% st_transform('ESRI:102658')

ggplot() +
  geom_sf(data=miami.base, fill="black") +
  geom_sf(data=worship, colour="red", size=.75)

MiamiDFKnown <-
  MiamiDFKnown %>% 
  mutate(
    worship_nn1 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(worship), 1),
    worship_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(worship), 2), 
    worship_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(worship), 10)) 

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, AllSchool, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2, worship_nn3))
summ(reg)
summary(reg)

## Parking

parking <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature(key = 'amenity', value = c("parking", "parking_space")) %>%
  osmdata_sf()

parking <- 
  parking$osm_points %>%
  .[miami.base,]

parking <- parking %>% st_transform('ESRI:102658')

ggplot() +
  geom_sf(data=miami.base, fill="black") +
  geom_sf(data=parking, colour="red", size=.75)

MiamiDFKnown <-
  MiamiDFKnown %>% 
  mutate(
    parking_nn1 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(parking), 1),
    parking_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(parking), 5), 
    parking_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(parking), 10)) 

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, AllSchool, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2,
                          worship_nn3, parking_nn1, parking_nn2, parking_nn3))
summ(reg)
summary(reg)

## Work centers
                  
landuse <- st_read("https://opendata.arcgis.com/datasets/244e956692d442c3beaa8a89259e3bd9_0.geojson")
landuse <- st_transform(landuse,'ESRI:102658')

office <- filter(landuse, DESCR == "Office Building.")

MiamiDFKnown <-
  MiamiDFKnown %>% 
  mutate(
    office_nn1 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(st_centroid(office)), 1),
    office_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(st_centroid(office)), 5), 
    office_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(st_centroid(office)), 10)) 

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, AllSchool, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2, 
                          worship_nn3, parking_nn1, parking_nn2, parking_nn3, office_nn1, office_nn2, office_nn3, CoastDist))
summ(reg)
summary(reg)


