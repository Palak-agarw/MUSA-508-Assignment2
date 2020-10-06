# Load Libraries
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
library(fastDummies)

options(scipen=999)
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

# Wrangle Miami Data
MiamiDF <- st_read("C:/Users/owner160829a/Documents/GitHub/MUSA-508-Assignment2/studentsData.geojson")

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

MiamiDF <- st_join(MiamiDF, Neighborhoods, join = st_intersects) 

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
  na.omit() %>% select("SalePrice", "XFs_Pool 8' res BETTER 3-8' dpth",
                       "XFs_Luxury Pool - Better","XFs_Luxury Pool - Good.",
                       "XFs_Pool - Wading - 2-4' depth",
                       "XFs_Pool COMM BETTER 3-6' dpth",
                       "XFs_Central A/C (Aprox 400 sqft/Ton)")

ggcorrplot(
  round(cor(dummyVars), 1), 
  p.mat = cor_pmat(dummyVars),
  colors = c("#25CB10", "white", "#FA7800"),
  type="lower",
  insig = "blank") +  
  labs(title = "Correlation across XF variables")   

# Regression
reg <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
             dplyr::select(SalePrice, Bed, Bath, Stories, YearBuilt,LivingSqFt))
summ(reg)
summary(reg)

# Adding Spatial Features
## Beaches (Looks weird)
Beaches<- st_read("https://opendata.arcgis.com/datasets/9e30807e3efd44f3b16ab8d3657249f2_0.geojson")

Beaches <- Beaches %>% 
  st_as_sf(coords = c("Lon", "Lat"), crs = 4326, agr = "constant") %>%
  st_transform('ESRI:102658')

ggplot() +
  geom_sf(data = Neighborhoods, fill = "grey40") +
  geom_sf(data = Beaches, aes(color="pink"), 
          show.legend = "point", size = 1) +
  mapTheme()

## Shoreline
Shoreline <- st_read("https://opendata.arcgis.com/datasets/58386199cc234518822e5f34f65eb713_0.geojson")
Shoreline<- st_transform(Shoreline,'ESRI:102658')

ggplot() +
  geom_sf(data = Neighborhoods, fill = "grey40") +
  geom_sf(data = Shoreline, color = "blue")

MiamiDFKnown <- mutate(MiamiDFKnown, distancetoshore=(st_distance(MiamiDFKnown, Shoreline)))

## Something weird going on with this regression
regshoreline <- lm(SalePrice ~ ., data = st_drop_geometry(MiamiDFKnown) %>% 
            dplyr::select(SalePrice, distancetoshore))
summ(regshoreline)
summary(regshoreline)

# School Attendance Areas
## Doesn't include Miami Beach
## Something is wrong with projection
ElementarySchool <- st_read("https://opendata.arcgis.com/datasets/19f5d8dcd9714e6fbd9043ac7a50c6f6_0.geojson")
ElementarySchool<- st_transform(ElementarySchool,'ESRI:102658') %>%
  select(-FID,-ID,-ZIPCODE,-PHONE,-REGION,-ID2,-FLAG,-CREATEDBY,-CREATEDDATE,-MODIFIEDBY,-MODIFIEDDATE)
ElementarySchool<-filter(ElementarySchool,CITY==c("Miami","MiamiBeach"))

ggplot() +
  geom_sf(data = ElementarySchool, fill = "grey40") +
  geom_sf(data = MiamiDF)+
  geom_sf(data = Neighborhoods)
  
MiamiDFSchool<-st_join(MiamiDF, ElementarySchool, join = st_intersects)
