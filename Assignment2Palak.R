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

# Data Wrangling

MiamiDF <- st_read("studentsData.geojson")

MiamiTrainingDF <- MiamiDF[!(MiamiDF$SalePrice==0),]

MiamiDF <- st_transform(MiamiTrainingDF,'ESRI:102658') #%>% 
#select(-saleQual,-GPAR,-Land.Use,-Owner1,-Owner2,-saleDate,-Year)

MiamiTrainingDF <- MiamiTrainingDF %>% st_transform('ESRI:102658')

Neighborhoods <- st_read("https://opendata.arcgis.com/datasets/2f54a0cbd67046f2bd100fb735176e6c_0.geojson")%>%
  st_transform('ESRI:102658')

metroStops <- st_read("https://opendata.arcgis.com/datasets/ee3e2c45427e4c85b751d8ad57dd7b16_0.geojson") 
metroStops <- metroStops %>% st_transform('ESRI:102658')


# Plot of the metro stops
ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data=MiamiTrainingDF)+
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
  st_intersection(buffer, MiamiTrainingDF) %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Clip")

ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data = clip) 

# Do a spatial selection to see which tracts touch the buffer
selection <- 
  MiamiTrainingDF[buffer,] %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Spatial Selection")

ggplot() +
  geom_sf(data=Neighborhoods)+
  geom_sf(data = selection) +
  geom_sf(data = buffer, fill = "transparent", color = "red")

selectCentroids <-
  st_centroid(MiamiTrainingDF)[buffer,] %>%
  st_drop_geometry() %>%
  left_join(dplyr::select(MiamiTrainingDF, Folio)) %>%
  st_sf() %>%
  dplyr::select(Folio) %>%
  mutate(Selection_Type = "Select by Centroids")

ggplot() +
  geom_sf(data = selectCentroids) +
  geom_sf(data = buffer, fill = "transparent", color = "red")+
  labs(title="Selected Tracts with centroids in the buffer zone", 
       subtitle="Bay Area, CA", 
       caption="Figure 2.2") +
  theme(plot.title = element_text(size=22))

MiamiTrainingDF.t <- 
  rbind(
    st_centroid(MiamiTrainingDF)[buffer,] %>%
      st_drop_geometry() %>%
      left_join(MiamiTrainingDF) %>%
      st_sf() %>%
      mutate(TOD = 1),
    st_centroid(MiamiTrainingDF)[buffer, op = st_disjoint] %>%
      st_drop_geometry() %>%
      left_join(MiamiTrainingDF) %>%
      st_sf() %>%
      mutate(TOD = 0))

MiamiTrainingDF <- MiamiTrainingDF %>%
  transform(TOD = ifelse(selection$Folio == MiamiTrainingDF$Folio,1,0))

# Correlation Matrix

numericVars <- 
  select_if(st_drop_geometry(MiamiTrainingDF.t), is.numeric) %>% na.omit() %>%
  select(SalePrice,Bed, Bath, Stories, Units, YearBuilt, LivingSqFt, TOD)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 

# Regression
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiTrainingDF) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt))
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

MiamiTrainingDF <- MiamiTrainingDF %>%
  mutate(geomtry = st_centroid(MiamiTrainingDF))

## Nearest Neighbor Feature
st_c <- st_coordinates

MiamiTrainingDF <-
  MiamiTrainingDF %>% 
  mutate(
    bar_nn1 = nn_function(st_c(st_centroid(MiamiTrainingDF)), st_c(bars), 1),
    bar_nn2 = nn_function(st_c(st_centroid(MiamiTrainingDF)), st_c(bars), 2), 
    bar_nn3 = nn_function(st_c(st_centroid(MiamiTrainingDF)), st_c(bars), 3), 
    bar_nn4 = nn_function(st_c(st_centroid(MiamiTrainingDF)), st_c(bars), 4), 
    bar_nn5 = nn_function(st_c(st_centroid(MiamiTrainingDF)), st_c(bars), 5)) 

## bar cor
st_drop_geometry(MiamiTrainingDF) %>% 
  mutate(Age = 2015 - YearBuilt) %>%
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
  select(SalePrice,Bed, Bath, Stories, Units, YearBuilt, LivingSqFt, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5)

ggcorrplot(
  round(cor(numericVars), 1), 
  p.mat = cor_pmat(numericVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across numeric variables") 

# Regression
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiTrainingDF) %>% 
            dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt, bar_nn1, bar_nn2, bar_nn3, bar_nn4, bar_nn5))
summ(reg)
summary(reg)