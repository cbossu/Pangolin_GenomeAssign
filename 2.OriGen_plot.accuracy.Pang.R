library(tidyverse)
library(geosphere)
library(beeswarm)
library(rworldmap)
library(WGCNA)
library(reshape2)
library(xlsx)

t<-read_delim("data/OriGen_input/Pang.geo_ref.combined.clip.ovlp.predictLatLong.txt",delim="\t")# %>% filter(Sample_Name %in% ind111$Sample)

ind<-read_delim("data/OriGen_input/Pangolin.111breeding_b.4Origen.loc",delim="\t") %>% rename(Sample_Name=Sample)
cv<-read.xlsx("data/Pangolin_sample_trackerv17.08182021.xlsx",sheetIndex=1) %>%  dplyr::select(Sample_name,Region_K5) %>% rename(Sample_Name=Sample_name)

compare<-t %>% 
  left_join(ind) %>% 
  left_join(cv)
compare$Sample_Name
compare$predLong
compare %>% dplyr::select(Region_K5) %>% distinct

mat <- distm(compare[,c('predLong','predLat')], compare[,c('Longitude','Latitude')], fun=distHaversine)
colnames(mat) <- compare$Sample_Name
rownames(mat) <- compare$Sample_Name

km1<-melt(mat) %>% filter(Var1==Var2) %>% rename(Sample_Name=Var1) %>% left_join(compare)

par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,0,0),tck=0.01)
par(mfrow=c(1,1))
pdf('results/OriGen_results/Pangolin.nofilt.OriGen_clipped.ovlp.grey.km2.RubCer.pdf')
hist(km1$value/1000,main="",xlab="Distance Error (km)",breaks=20,col="gray40",border="white")
dev.off()

hist(compare$highprob,main="",xlab="High prob histogram",breaks=20,col="gray40",border="white")

compare_filt<- compare %>% filter(highprob>0.05)

mat <- distm(compare_filt[,c('predLong','predLat')], compare_filt[,c('Longitude','Latitude')], fun=distHaversine)
colnames(mat) <- compare_filt$Sample_Name
rownames(mat) <- compare_filt$Sample_Name

km1<-melt(mat) %>% filter(Var1==Var2) %>% rename(Sample_Name=Var1) %>% left_join(compare)

rub<-read_delim("results/rubias_results/Pangolin.self_assessment.rubias.InferGroup.txt",delim="\t") %>% rename(Sample_Name=indiv) %>% dplyr::select(Sample_Name,repunit,inferred_repunit,collection)

km1 %>% left_join(rub)  %>% dplyr::select(Sample_Name,collection,repunit,inferred_repunit,Region_K5,Latitude,Longitude,everything()) %>% mutate(error_km=value/1000) %>% dplyr::select(-Var2,-value) %>% write.table("results/OriGen_results/Pangolin.Origen.self_assessment.rubias.InferGroup.Error2.txt",row.names=F,quote=F,sep="\t")


err_origen<-read_delim("results/rubias_results/Pangolin.self_assessment.rubias.InferGroup.Error2.txt",delim="\t") %>% 
  filter(Sample_Name %in% ind111$Sample) %>% 
  dplyr::select(-Region_K5) %>% 
  filter(error_km>500) %>% 
  #filter(repunit!=inferred_repunit) %>% 
  #filter(Sample_Name=="CA_ECO88094LE")


par(mfrow=c(1,2),mar=c(3,3,1,1),mgp=c(1.5,0,0),tck=0.01)
par(mfrow=c(1,1))
pdf('Pangolin.filt_highprob.OriGen_clipped.ovlp.grey.km2.pdf')
hist(km1$value/1000,main="",xlab="Distance Error (km)",breaks=20,col="gray40",border="white")
dev.off()


my_cols <- paste(c('#FFCC00','#CC0000','#33FFFF','#9933CC','#009933'),"80",sep="")
newmap <- getMap(resolution="low")

plot(1:10,1:10,type="n",xlim=c(-15,40),ylim=c(-15,15),
     xlab="Longitude",ylab="Latitude")
plot(newmap,lwd=1.5,border="gray30",add=T,col="gray90")
points(compare$predLong,compare$predLat,col=my_cols,pch=16,cex=0.9)

breeding <- st_read("~/Dropbox/BGP/Pangolin/Phataginus_tricuspis/data_0.shp")

library(rgdal)
my_spdf<-st_read("~/Dropbox/BGP/Pangolin/meta_data/world_shape_file/TM_WORLD_BORDERS_SIMPL-0.3.shp")
# Keep only data concerning Africa
africa <- my_spdf %>% filter(REGION==2)

africa_geo = st_transform(africa, 4326)
extent(africa_geo) 

domain <- c(
  xmin = -220, 
  xmax = -20,
  ymin = 15,
  ymax = 95
)


coastlines <- st_read("~/Dropbox/BGP/genoscape_maps/shapefiles/ne_shapefiles/ne_10m_coastline.shp")
coast_cropped <- st_crop(coastlines, africa_geo)
countries_cropped <-  st_read("~/Dropbox/BGP/genoscape_maps/shapefiles/ne_shapefiles/ne_10m_admin_0_boundary_lines_land.shp") %>%
  st_crop(africa_geo)
km_high<-km1 %>% filter(value/1000 > 1000)
library(ggspatial)
rectangled <- ggplot() +
  ggspatial::layer_spatial(coast_cropped) +
  geom_sf(data = countries_cropped, fill = NA, size = 0.15) +
  #ggspatial::layer_spatial(genoscape_rgba) +
  geom_spatial_point(data = km1, mapping = aes(x = predLong, y = predLat, col=Region_K5),size=1) +
  geom_spatial_text(data = km_high, mapping = aes(x = predLong, y = predLat,label=Sample_Name),size=2) +
  geom_sf(data=breeding,fill=NA,size=.5)  +
  scale_color_manual(values = my_cols) +
  theme_bw() + theme(legend.key.size = unit(0.1, "cm"))+
  coord_sf(xlim = c(-16.7559, 35.5713),
           ylim = c(-15.8594, 15.3875),
           expand = FALSE) +
  xlab("Longitude")+ylab("Latitude")
rectangled

ggsave("Pangolin.OriGen_accuracy.nofilt.outline.clip.pdf")

km1<-read_delim("results/OriGen_results/Pangolin.Origen.self_assessment.rubias.InferGroup.Error2.txt")
mean(km1$error_km,na.rm=T)
median(km1$error_km,na.rm=T)
km2<-read_delim("results/OriGen_results/Pangolin.Origen.self_assessment.rubias.InferGroup.Error2.txt") %>% filter(error_km<500)
mean(km2$error_km)
median(km2$error_km)
