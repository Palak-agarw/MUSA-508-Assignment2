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

MiamiDF <- st_as_sf(MiamiDF)

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

# Create Categorical Variables
MiamiDF$YearCat <- cut(MiamiDF$YearBuilt, c(1900,1909,1919,1929,1939,1949,
                                            1959,1969,1979,1989,1999,2009,2019))

MiamiDF$BedCat <- cut(MiamiDF$Bed,breaks=c(0:8,Inf), right=FALSE, labels=c(0:7,"8+"))

# Join Neighborhood Data
Neighborhoods <- st_read("https://opendata.arcgis.com/datasets/2f54a0cbd67046f2bd100fb735176e6c_0.geojson")

Neighborhoods <- st_transform(Neighborhoods,'ESRI:102658')

Municipality <- st_read('https://opendata.arcgis.com/datasets/bd523e71861749959a7f12c9d0388d1c_0.geojson')

Municipality <- st_transform(Municipality,'ESRI:102658')

Municipality <- filter(Municipality, NAME == "MIAMI BEACH")

Neighborhoods <- Neighborhoods %>%
  rename( Neighbourhood_name = LABEL)

Neighborhoods <- Neighborhoods %>% dplyr::select(-FID,-Shape_STAr,-Shape_STLe,-Shape__Area,-Shape__Length)

Municipality <- Municipality %>%
  rename( Neighbourhood_name = NAME) 
  
Municipality <- Municipality %>% dplyr::select(-FID,-MUNICUID,-MUNICUID,-FIPSCODE,-CREATEDBY, -CREATEDDATE, -MODIFIEDBY, 
                                        -MODIFIEDDATE, -SHAPE_Area, -SHAPE_Length, -MUNICID)


Municipality$Neighbourhood_name <- make.names(Municipality$Neighbourhood_name, unique=TRUE)

Neighborhoods_combine <- rbind(Neighborhoods, Municipality)

MiamiDF <- st_as_sf(MiamiDF)

MiamiDF <-   st_centroid(MiamiDF)

MiamiDF <- st_join(MiamiDF, Neighborhoods_combine, join = st_intersects) 


## Census data

census_api_key("aea3dee2d96acb5101e94f3dcfa1b575f73d093a", overwrite = TRUE)

Miamitracts <-  
  get_acs(geography = "tract", variables = c("B25026_001E","B02001_002E","B19013_001E",
                                             "B25058_001E","B02001_003E","B02001_004E",
                                             "B02001_005E","B02001_006E","B03001_003E"), 
          year=2017, state=12, county=086, geometry=T) %>% 
  st_transform('ESRI:102658')

Miamitracts <- 
  Miamitracts %>%
  dplyr::select( -NAME, -moe) %>%
  spread(variable, estimate) %>%
  dplyr::select(-geometry) %>%
  rename(TotalPop = B25026_001, 
         Whites = B02001_002,
         Blacks = B02001_003,
         AmInd = B02001_004,
         Asian = B02001_005,
         Hawaiian = B02001_006,
         Hispanic = B03001_003,
         MedHHInc = B19013_001, 
         MedRent = B25058_001
  )

Miamitracts <-
  Miamitracts %>%
  mutate(pctWhite = ifelse(TotalPop > 0, Whites / TotalPop * 100, 0),
         pctHispanic = ifelse(TotalPop > 0, Hispanic / TotalPop * 100,0))


ggplot() + 
  geom_sf(data=Miamitracts)

MiamiDF <- st_join(MiamiDF, Miamitracts, join = st_intersects) 

# Remove Challenge Houses

MiamiDFKnown <- MiamiDF[!(MiamiDF$SalePrice==0),]

# Map House Sales and Neighborhoods

ggplot() +
  geom_sf(data=Neighborhoods)+
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

# Correlation Matrices
numericVars <- 
  select_if(st_drop_geometry(MiamiDFKnown), is.numeric) %>% na.omit() %>%
  select(SalePrice, Land, Bldg, Total, Assessed, City.Taxable, AdjustedSqFt,
         LotSize, Bed, Bath, Stories, Units, YearBuilt, EffectiveYearBuilt,
         LivingSqFt, ActualSqFt)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 

dummyVars <- select_if(st_drop_geometry(MiamiDFKnown), is.numeric) %>% 
  na.omit() %>% select("SalePrice", "3to8ftPool",
                       "3to6ftPool","LuxuryPool",
                       "8ftres3to8ftPool",
                       "Whirpool","2to4ftPool",
                       "XFs_Central A/C (Aprox 400 sqft/Ton)")

ggcorrplot(
  round(cor(dummyVars), 1), 
  p.mat = cor_pmat(dummyVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across XF variables")   

# Regression Adj R-squared=..7276
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt))
summ(reg)
summary(reg)

## Regression with Neighborhoods Adj R-squared=.7474
reg2 <-lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt,Neighbourhood_name))
summ(reg2)
summary(reg2)

# Creating Year Built Categories A rsquared=.7539
## Need to move this code before split

reg3 <-lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt,Neighbourhood_name))
summ(reg3)
summary(reg3)


#Transit Data 

metroStops <- st_read("https://opendata.arcgis.com/datasets/ee3e2c45427e4c85b751d8ad57dd7b16_0.geojson") 
metroStops <- metroStops %>% st_transform('ESRI:102658')


## Buffer method - not preferred
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
    Distance >= 0.5 & Distance <= 0.75  ~ "Three Quater Mile",
    Distance >= 1 & Distance < 2   ~ "More than one Mile",
    Distance >= 2 & Distance < 3   ~ "More than two Mile",
    Distance >= 3    ~ "More than three Mile"))

# Regression R-sq = .7695

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, Neighbourhood_name, NewDistance.cat))
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

# Regression R-sq = .7732
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, Neighbourhood_name,
                          NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5))
summ(reg)
summary(reg)

## Shoreline

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

## Regression R-sq = .7738
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
                     dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, Neighbourhood_name, 
                                   NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5, CoastDist))
summ(reg)
summary(reg)

# School Attendance Areas
## Elementary Schools
## Doesn't include Miami Beach
## Something is wrong with projection
ElementarySchool <- st_read("https://opendata.arcgis.com/datasets/19f5d8dcd9714e6fbd9043ac7a50c6f6_0.geojson")

ElementarySchool<- st_transform(ElementarySchool,'ESRI:102658') %>%
  dplyr::select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-FLAG,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

ElementarySchool<-filter(ElementarySchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = ElementarySchool, fill = "grey40") +
  geom_sf(data = MiamiDF)+
  geom_sf(data = Neighborhoods)

## Middle Schools
MiddleSchool <- st_read("https://opendata.arcgis.com/datasets/dd2719ff6105463187197165a9c8dd5c_0.geojson")

MiddleSchool<- st_transform(MiddleSchool,'ESRI:102658') %>%
  dplyr::select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

MiddleSchool<-filter(MiddleSchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = MiddleSchool, fill = "grey40") +
  geom_sf(data = MiamiDF)


## High Schools
HighSchool <- st_read("https://opendata.arcgis.com/datasets/9004dbf5f7f645d493bfb6b875a43dc1_0.geojson")

HighSchool<- st_transform(HighSchool,'ESRI:102658') %>%
  dplyr::select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)

HighSchool<-filter(HighSchool, CITY == "Miami"| CITY == "MiamiBeach")

ggplot() +
  geom_sf(data = HighSchool, fill = "grey40") +
  geom_sf(data = MiamiDF) +
  geom_sf(data = Neighborhoods)

ElementarySchool['Elementary'] <- "Elementary"
MiddleSchool['Middle'] <- "Middle"
HighSchool['High']<-"High"

ElementarySchool <- ElementarySchool %>% 
  dplyr::select(-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area)

MiddleSchool <- MiddleSchool %>% 
  dplyr::select(-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area)

HighSchool <- HighSchool %>% 
  dplyr::select(-ADDRESS,-CITY,-GRADES,-DISPLAYNAME,-SHAPE_Length,-SHAPE_Area)

MiamiDFKnown <- st_join(MiamiDFKnown, ElementarySchool, join = st_intersects) 
MiamiDFKnown <- st_join(MiamiDFKnown, MiddleSchool, join = st_intersects) 
MiamiDFKnown <- st_join(MiamiDFKnown, HighSchool, join = st_intersects) 

MiamiDFKnown$NAME <- ifelse(is.na(MiamiDFKnown$NAME), "Other", MiamiDFKnown$NAME)
MiamiDFKnown$NAME.x <- ifelse(is.na(MiamiDFKnown$NAME.x), "Other", MiamiDFKnown$NAME.x)
MiamiDFKnown$NAME.y <- ifelse(is.na(MiamiDFKnown$NAME.y), "Other", MiamiDFKnown$NAME.y)

# Regression R-sq = .7738

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
                     dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt,Neighbourhood_name, NewDistance.cat, bar_nn1, 
                                   bar_nn2, bar_nn3, bar_nn4,bar_nn5, CoastDist, NAME.x, NAME, NAME.y))
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

# Regression R-sq = .7828

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt,Neighbourhood_name, 
                          NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5, CoastDist, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5 ))
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

# Regression R-sq = .7827
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2, 
                          worship_nn3, Neighbourhood_name))
summ(reg)
summary(reg)

## Parking

parking <- opq(bbox = c(xmin, ymin, xmax, ymax)) %>% 
  add_osm_feature(key = 'amenity', value = c("parking", "parking_space", "parking_entrance")) %>%
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
    parking_nn2 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(parking), 2), 
    parking_nn3 = nn_function(st_c(st_centroid(MiamiDFKnown)), st_c(parking), 5)) 

# Regression R-sq = .7826
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist, Neighbourhood_name, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2,
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

# regression R-sq = .783

reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>%
            dplyr::select(SalePrice, Bed, Bath, Stories, YearCat,LivingSqFt, NewDistance.cat, bar_nn1, bar_nn2, bar_nn3, bar_nn4, 
                          bar_nn5, CoastDist,Neighbourhood_name, park_nn1, park_nn2, park_nn3, park_nn4, park_nn5, worship_nn1, worship_nn2, 
                          worship_nn3, parking_nn1, parking_nn2, parking_nn3, office_nn1, office_nn2, office_nn3, CoastDist))
summ(reg)
summary(reg)

## Lag price

k_nearest_neighbors = 5
MiamiDFKnownDistinct <- distinct(MiamiDFKnown,geometry,.keep_all=TRUE)
coords <- st_coordinates(MiamiDFKnownDistinct)
neighborList <- knn2nb(knearneigh(coords, k_nearest_neighbors))
spatialWeights <- nb2listw(neighborList, style="W")
MiamiDFKnownDistinct$lagPrice <- lag.listw(spatialWeights, MiamiDFKnownDistinct$SalePrice)

MiamiDFKnownDistinct <- dplyr::select(MiamiDFKnownDistinct,geometry,lagPrice)
MiamiDFKnown <- st_join(MiamiDFKnown,MiamiDFKnownDistinct,left=TRUE)

# External Model Validation
## set random seed
set.seed(121491)

# get index for training sample
inTrain <- caret::createDataPartition( 
  y = paste(MiamiDFKnown$Neighbourhood_name,MiamiDFKnown$BedCat, MiamiDFKnown$YearCat), 
  p = .60, list = FALSE)

# split data into training and test, before comma is row, after comma is column
Miami.training <- MiamiDFKnown[inTrain,] 
Miami.test     <- MiamiDFKnown[-inTrain,]  

# Regression  0.7862, 0.7927, 0.7882, 0.7835
reg5 <- lm(SalePrice ~ ., data = st_drop_geometry(Miami.training) %>% 
             dplyr::select(SalePrice,BedCat, Bath, Stories, YearCat,LivingSqFt,Patio,Fence,`3to8ftPool`,Whirpool,
                           `3to6ftPool`,`LuxuryPool`, NewDistance.cat, bar_nn2, bar_nn3, bar_nn4, Neighbourhood_name,
                           CoastDist, park_nn4, park_nn5, CoastDist, lagPrice.y, MedHHInc, LotSize, ActualSqFt))


summary(reg5)

## predicting on new data
reg5_predict <- predict(reg5, newdata = Miami.test)

nhood_sum <- Miami.test %>% 
  group_by(Neighbourhood_name) %>%
  summarize(meanPrice = mean(SalePrice, na.rm = T),
            meanPrediction = mean(SalePrice.Predict, na.rm = T),
            meanMAE = mean(SalePrice.AbsError, na.rm = T))

nhood_sum %>% 
  st_drop_geometry %>%
  arrange(desc(meanMAE)) %>% 
  kable() %>% kable_styling()

Miami.test <-
  Miami.test %>%
  mutate(SalePrice.Predict = predict(reg5, Miami.test),
         SalePrice.Error = SalePrice.Predict - SalePrice,
         SalePrice.AbsError = abs(SalePrice.Predict - SalePrice),
         SalePrice.APE = (abs(SalePrice.Predict - SalePrice)) / SalePrice.Predict)

mean(Miami.test$SalePrice.AbsError, na.rm = T)
mean(Miami.test$SalePrice.APE, na.rm = T)

### Try to figure out how to make a histogram of AbsError, see section 3.4.1 ## add different x value
hist(Miami.test$`SalePrice.AbsError`,xlab="Sales Price Absolute Error",
     breaks=50, col="purple")

# Measure Generalizability
## Mean Square Error train and test
rmse.train <- caret::MAE(predict(reg5), Miami.training$SalePrice)
rmse.test  <- caret::MAE(reg5_predict, Miami.test$SalePrice)

cat("Train MAE: ", as.integer(rmse.train), " \n","Test MAE: ", as.integer(rmse.test))

# Plotting accuracy metrics
preds.train <- data.frame(pred   = predict(reg5),
                          actual = Miami.training$SalePrice,
                          source = "training data")
preds.test  <- data.frame(pred   = reg5_predict,
                          actual = Miami.test$SalePrice,
                          source = "testing data")
preds <- rbind(preds.train, preds.test)

ggplot(preds, aes(x = pred, y = actual, color = source)) +
  geom_point() +
  geom_smooth(method = "lm", color = "green") +
  geom_abline(color = "orange") +
  coord_equal() +
  theme_bw() +
  facet_wrap(~source, ncol = 2) +
  labs(title = "Comparing predictions to actual values",
       x = "Predicted Value",
       y = "Actual Value") +
  theme(
    legend.position = "none"
  )

# Cross Validation
fitControl <- trainControl(method = "cv", 
                           number = 10,
                           # savePredictions differs from book
                           savePredictions = TRUE)

set.seed(717)

# for k-folds CV
reg.cv <- 
  train(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
          dplyr::select(SalePrice,BedCat, Bath, Stories, YearCat,LivingSqFt,Patio,Fence,`3to8ftPool`,Whirpool,
                        `3to6ftPool`,`LuxuryPool`, NewDistance.cat, bar_nn2, bar_nn3, bar_nn4, Neighbourhood_name,
                        CoastDist, park_nn4, park_nn5, office_nn2, office_nn3, CoastDist, lagPrice), 
        method = "lm", 
        trControl = fitControl, 
        na.action = na.pass)

reg.cv

reg.cv$resample

reg.cv$resample %>% 
  pivot_longer(-Resample) %>% 
  mutate(name = as.factor(name)) %>% 
  ggplot(., aes(x = name, y = value, color = name)) +
  geom_jitter(width = 0.1) +
  facet_wrap(~name, ncol = 3, scales = "free") +
  theme_bw() +
  theme(
    legend.position = "none"
  )

# extract predictions from CV object
cv_preds <- reg.cv$pred
# compare number of observations between data sets
nrow(MiamiDFKnown)
nrow(cv_preds)

## Create dataset with "out of fold" predictions and original data

MiamiDFKnown <- MiamiDFKnown %>%
  mutate(pred = predict(reg.cv, .)) %>%  #get a prediction for each row
  mutate(error = pred - SalePrice)    #calculate errors

map_preds <- MiamiDFKnown %>% 
  rowid_to_column(var = "rowIndex") %>% 
  left_join(cv_preds, by = "rowIndex") %>% 
  mutate(SalePrice.AbsError = abs(pred - SalePrice)) %>% 
  cbind(st_coordinates(.))

# weird CRS fix to boston.sf
st_crs(map_preds) <- st_crs(Neighborhoods)

# plot errors on a map
ggplot() +
  geom_sf(data = Neighborhoods_combine, fill = "grey40") +
  geom_sf(data = MiamiDFKnown, aes(colour = q5(error)),
          show.legend = "point", size = 1) +
  scale_colour_manual(values = palette5,
                      labels=qBr(MiamiDFKnown,"error"),
                      name="Quintile\nBreaks") +
  labs(title="Absolute sale price errors on the OOF set",
       subtitle = "OOF = 'Out Of Fold'") +
  mapTheme()

map_preds_sum <- MiamiDFKnown %>% 
  group_by(Neighbourhood_name) %>% 
  summarise(meanMAE = mean(error))

ggplot() +
  geom_sf(data = Neighborhoods_combine %>% 
            left_join(st_drop_geometry(map_preds_sum), by = "Neighbourhood_name"),
          aes(fill = q5(meanMAE))) +
  scale_fill_manual(values = palette5,
                    labels=qBr(nhood_sum,"meanMAE"),
                    name="Quintile\nBreaks") +
  mapTheme() +
  labs(title="Absolute sale price errors on the OOF set by Neighborhood")
