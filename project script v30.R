library(factoextra)
library(readxl)
library(writexl)
library(vegan)
library(geosphere)
library(terra)
library(ggfortify)
library(RNetCDF) 
library(RStoolbox)
library(ncmeta)
library(stars)
library(gdalcubes)
library(ncdf4)
library(sf)
library(mapsf)
library(dplyr)
library(mi)
library(scales)
library(plotly)
library(psych)

library(ade4)
library(ggplot2)
library(dichromat)
library(tidyverse)


setwd("D://2. Hyperspectral redbluff Analysis/")


###**Step 1: Requirements (orthorectified geotiff and .nc file from EMIT)**
##The only variable thats especially important from the .nc file is Wvl for plotting later

nc_files <- NDVI_files <- list.files("./EMIT files/", pattern = ".nc", full.names = T)
nc <- nc_files[[3]] # 1 = 04/15/23, 2= 10/07/24, 10/15/24, all should work

#open ncdf file to extract variables
nc_file <- ncdf4::nc_open(nc)

print(ncdf4::ncatt_get(nc_file, 0))
print(nc_file)
# Extract dimensions and the variable name (adjust as needed)
data_var <- "reflectance" # Replace with actual variable name
wvl <- ncvar_get(nc_file, "wavelengths") # wvls
data_cube <- ncvar_get(nc_file, data_var)  # 3D array: [lon, lat, layer]
#Spatial Variables
Crs <- terra::crs("epsg:4326")
bands <- c(1:length(wvl))
bands_wvl <- paste0(bands, "_", wvl)
#.nc file is no longer useful for this code 
ncdf4::nc_close(nc_file)
rm(data_cube)
gc()


bands
wvl<- substr(wvl, 1, 6)
bands_wvl


###**Step 2: Import geotiff, project correctly before PCA**

#Import geotiff !!!!rename
tif_files <- NDVI_files <- list.files("./Inputs/", pattern = "ortho_geotiff.tif", full.names = T)
extracted_name <- basename(sub("_.*", "", tif_files))
output_folder <- "./Outputs/" 
refl_stack <- c()
for (x in seq_along(tif_files)) {
  fp <- terra::rast(tif_files[x]) #rasterizes tifs in input folder
  id <- extracted_name[x] #used for naming each layer
  rb.shp <- terra::vect("./Inputs/site shape files/RedBluffBoundary.shp") #cropping to redbluff exe
  rb.shp <- terra::project(rb.shp, crs(fp))
  fp_crop <- terra::crop(fp, rb.shp)
  refl_tif <- terra::mask(fp_crop, rb.shp)
  names(refl_tif) <- bands_wvl
  
  name <- paste0("refl_tif_", id, "full_spec")
 # NIR_name <- paste0("refl_tif_", id, "NIR") unused for now
 # RGB_name <- paste0("refl_tif_", id, "RGB")
  
  assign(name, refl_tif)
  
  refl_tif_sub <- refl_tif[[1]]
  
  if(is.null(refl_stack)) {refl_stack <- name} else { refl_stack <- c(refl_stack, name)}
 # if(is.null(sub_refl_stack)) {sub_refl_stack <- refl_tif_sub } else { sub_refl_stack <- c(sub_refl_stack, refl_tif_sub)}
  #creating a RGB image
  red <- fp[[40]]
  green <- fp[[21]]
  blue <- fp[[13]]
  #nir_stack<- fp[[]]
  #nir <- app(nir_stack, fun=mean, na.rm=TRUE)
  #swir_stack <- fp[[]]
  #nir <- app(nir_stack, fun=mean, na.rm=TRUE)
  
   
  rgb_stack <- c(blue, green, red)
  terra::plotRGB(rgb_stack, r=1, g=2, b=3, stretch = "lin")
  
  
  writeRaster(rgb_stack, filename = paste0(output_folder, "RGB_", id, ".tif"), filetype= "GTiff", overwrite = TRUE)
 
}
terra::plot(refl_tif[[26]])
rm(refl_tif)#removes extra output

###**Step 2.1: import snow raster**

snow_persistance <- rast("./inputs/rb_ndsi_snowper.tif") #snow persistance raster
snow_persistance <- project(snow_persistance, rb.shp)
snow_crop <- terra::crop(snow_persistance, rb.shp)
snow_persistance <- terra::mask(snow_crop, rb.shp) #snow cropped and masked
plot(snow_persistance, main="snow_persistance")

###**Step 3: PCA time**
anthro_mask <- terra::vect("./Inputs/site shape files/rb_anthro_mask.shp")
anthro_mask <- project(anthro_mask, crs(snow_persistance))

for (x in seq_along(refl_stack)){
  raster_object <- get(refl_stack[[x]])
  id <- extracted_name[x]
  raster_object_masked <- terra::mask(raster_object, anthro_mask, inverse=T)
  ras_df <- as.data.frame(raster_object_masked, xy = TRUE, na.rm = TRUE)
  colnames(ras_df) <- c("x", "y", wvl)
  ras_dfw <- ras_df[, sapply(1:ncol(ras_df), function(i) length(unique(ras_df[,i]))) != 1] #filters columns with no variation
  ras_xy <- ras_df %>% select(1:2) #grabs coordinates for rasterization later
  ras_dfw <- ras_dfw %>% select(3:ncol(.)) #removes xy before pca
  pcs <- principal(ras_dfw, rotate = "none", nfactors = 10, scores = TRUE) 
  #pcs <- principal(ras_dfw, rotate = "varimax", nfactors = 4, scores = TRUE) 
  loadings <- as.matrix(pcs$loadings)
  loadings_name <- paste0("loadings_", id)
  assign(loadings_name, loadings)
  
  pca_scores <- as.data.frame(scores(pcs))
  pca_scores <- bind_cols(ras_xy, pca_scores)
  
  #backtrans to raster
  pcs99 <- rast(pca_scores, type = 'xyz', crs = "EPSG:4326") 

  pcs_name <- paste0("pcs99_", id)
  assign(pcs_name, pcs99)
  
              # Save each plot
              # 1. PC1 vs PC2
              jpeg(filename = paste0(output_folder, "PCA_PC1_vs_PC2_", id, ".jpg"), width = 800, height = 600)
              plot(pcs99$PC1, pcs99$PC2, main = "PC1 vs PC2", xlab = "PC1", ylab = "PC2")
              dev.off()
              
              # 2. PC1 vs PC3
              jpeg(filename = paste0(output_folder, "PCA_PC1_vs_PC3_", id, ".jpg"), width = 800, height = 600)
              plot(pcs99$PC1, pcs99$PC3, main = "PC1 vs PC3", xlab = "PC1", ylab = "PC3")
              dev.off()
              
              # 3. PC1 vs PC4
              jpeg(filename = paste0(output_folder, "PCA_PC1_vs_PC4_", id, ".jpg"), width = 800, height = 600)
              plot(pcs99$PC1, pcs99$PC4, main = "PC1 vs PC4", xlab = "PC1", ylab = "PC4")
             dev.off()
              
              # 4. Histogram of PCs
              jpeg(filename = paste0(output_folder, "PCA_Histogram_", id, ".jpg"), width = 800, height = 600)
              hist(pcs99, main = "Histogram of PCA Components", xlab = "Component Value", ylab = "Frequency", col = "lightblue")
             dev.off()
              
              # 5. RGB Plot
              jpeg(filename = paste0(output_folder, "PCA_RGB_", id, ".jpg"), width = 800, height = 600)
              plotRGB(pcs99, r = 1, g = 2, b = 3, stretch = "hist", main = "RGB Representation of PCA")
              dev.off()
              writeRaster(pcs99, filename = paste0(output_folder, "PCA_RGB_", id, ".tif"), filetype= "GTiff", overwrite = TRUE)
}

plotRGB(pcs99_003,r = 1,g = 2,b = 3,stretch="hist") #trippy
plotRGB(pcs99_016,r = 1,g = 2,b = 3,stretch="hist") 
plotRGB(pcs99_018,r = 1,g = 2,b = 3,stretch="hist") 

##**Step: 3.2: interactive PCA plot**
# Extract the values for the first three principal components

pcs99 <- pcs99_018 #the below code depends on this line


pcs_data <- data.frame(
  PC1 = values(pcs99[[1]]),
  PC2 = values(pcs99[[2]]),
  PC3 = values(pcs99[[3]])
)

# Remove rows with NA values
pcs_data <- na.omit(pcs_data)

# Create a 3D scatter plot
fig <- plot_ly(
  data = pcs_data,
  x = ~PC1,
  y = ~PC2,
  z = ~PC3,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 2, color = pcs_data$PC1, colorscale = "Viridis")
)

# Display the plot
fig


###**Part 4: Import '23 field season data and X/Y file**
xy_23 <- read_excel("./Inputs/x,y and site names (edited spring25).xlsx", sheet = 1)
data_23 <- read_excel("./Inputs/2023 field season data (3-3-25).xlsx", sheet=4)

#checking xy data for NA's
xy_pts <- st_read('./Inputs/23 site xy shape files/Fieldseason_23_plotlocations.shp')

plot(xy_pts)
class(xy_pts)
str(xy_pts)
class(xy_pts)

#checking field data for NA's
dim(data_23)
data_missing <-data_23 %>% as.data.frame() %>% missing_data.frame()
image(data_missing)

#Cleaning Data and creating focal species subset

focal <- data_23 %>% dplyr::select(Plot_ID, #choose which columns will be used
                                   Longitude,
                                   Latitude,
                                   `S. densa crust Ab`,
                                   `BlueBunch Ab`,
                                   `Idaho Fescue Ab`,
                                   `Prairie June Grass Ab`,
                                   `Western Wheatgrass Ab`,
                                   `Blue Grama Ab`,
                                   `Needle-and-Thread Ab`,
                                   `Smooth Brome Ab`,
                                   `Cheatgrass Ab`,
                                   `Basin Wildrye Ab`,
                                   `%C (Forb)`,
                                   `%C (Woody)`,
                                   `%C (Cacti)`) %>% 
                      dplyr::rename(seda=`S. densa crust Ab`, #rename for simplicity
                                 pssp = `BlueBunch Ab`,
                                 feid = `Idaho Fescue Ab`,
                                 koma = `Prairie June Grass Ab`,
                                 pasm = `Western Wheatgrass Ab`,
                                 bogr = `Blue Grama Ab`,
                                 heco = `Needle-and-Thread Ab`,
                                 brin = `Smooth Brome Ab`,
                                 brte = `Cheatgrass Ab`,
                                 leci = `Basin Wildrye Ab`,
                                 forb = `%C (Forb)`,
                                 woody = `%C (Woody)`,
                                 cacti = `%C (Cacti)`) %>% 
                      mutate(across(seda:brte, ~ as.numeric(as.character(.))))

#diversity <- focal %>% mutate(across(seda:cacti, ~ (. / 100)*log(./100))) %>% #H' is plot wide Shannon-Weiner Diversity
#  mutate(across(seda:cacti, ~ tidyr::replace_na(., 0))) %>%
#  mutate("H'" = -rowSums(select(., seda:cacti), na.rm = TRUE))

                                  
topo <- data_23 %>% dplyr::select(Plot_ID, #choose which columns will be used
                                  Longitude,
                                  Latitude,
                                  Curvature_,
                                  Slope,
                                  Aspect,
                                  Elevation,
                                  TPI_90m,
                                  `%C (bare ground)`) %>%
                           rename(curve = 'Curvature_',
                                  aspect = 'Aspect',
                                  tpi = 'TPI_90m',
                                  bareground = `%C (bare ground)`) %>%
                           mutate(northness = cos(aspect*pi/180),
                                  eastness = sin(aspect*pi/180)) %>%
                          select(-aspect)
                           
#making site object for later
site <-c(rep(3,21),rep(5,21), rep(8,21), rep(9,21), rep(12,21), rep(18, 21), rep(22,21), rep(28,21), rep(40,21), rep(41, 21), rep(45, 21), rep(46,20), rep(49,21), rep(53,21), rep(54,21), rep(56,21), rep(59,21), rep(63, 21), rep(64,21), rep(65, 21), rep(66, 21))
site <- data.frame(site)
site <- site %>% dplyr::mutate(site = as.factor(site)) #should be treated as categorical value
dim(site)

###extracting snow persistence values to points
snow_points <- terra::extract(snow_persistance, xy_pts) 
snow_points <- cbind(xy_23,snow_points)

### extracting  pca to points
for(x in seq_along(extracted_name)){
  id <- extracted_name[[x]]
  pcs99_name <- paste0("pcs99_", id)
  pcs_99 <- get(pcs99_name)
  pc_points <- terra::extract(pcs_99, xy_pts) 
  xy_pc <- cbind(xy_23,pc_points)

  all_equal <- all(xy_pc$PlotID == xy_pc$ID)
  #checks that the two data sets line up after combining, breaks the loop if there's an issue
  if (all_equal) {
    print("The plot_id columns are equal in every row.")
  } else {
    print("The plot_id columns are NOT equal in every row.")
    break
  }
  
  
##**Part 5.5: subsetting and making dist matrices**##  
  
  ### Subset data (with unique id's for each data set)
  pcs_subset <- xy_pc %>% dplyr::select(PC1:PC10) 
  pcs_subset_name <- paste0("pcs_subset_", id)
  assign(pcs_subset_name, pcs_subset)
  
  snow_subset <- snow_points %>% dplyr::select(mean) %>%
              rename(snow_persistence = mean)
  
  focal_subset <- focal %>% dplyr::select(pssp:brte) #this sets veg vars (copy to below code if changed)
  
  
  topo_subset <- topo %>% select(curve:eastness) #this sets topo vars
  
  
  
  
  
  
  alldata <- cbind(site,pcs_subset, focal_subset, topo_subset, snow_subset) 
  alldata_name <- paste0("alldata_", id)
  assign(alldata_name, alldata)
  
  alldata_avg <- alldata %>% dplyr::group_by(site) %>%
    dplyr::summarize(across(everything(), mean, na.rm = TRUE)) #deprecated code still works tho
  alldata_avg_name <- paste0("alldata_avg_", id)
  assign(alldata_avg_name, alldata_avg)
  
  #re-subsetting new data set
  avg_pcs <- alldata_avg %>% dplyr::select(PC1:PC10) 
  avg_focal <- alldata_avg %>% dplyr::select(pssp:brte)
  avg_topo <- alldata_avg %>% dplyr::select(curve:snow_persistence) 
  avg_snow <- alldata_avg %>% dplyr::select(snow_persistence)
  
  #data is standardized (not between 0-1 though) (retry this )
  avg_pcs_scaled <- decostand(avg_pcs, method = "standardize")
  avg_focal_scaled <- decostand(avg_focal, method = "standardize")
  avg_topo_scaled <- decostand(avg_topo, method = "standardize")
  avg_snow_scaled <- decostand(avg_snow, method = "standardize")
    #giving final outputs unique ids
  avg_pcs_scaled_name <- paste0("avg_pcs_scaled_", id)
  assign(avg_pcs_scaled_name, avg_pcs_scaled)

}


###**Part 6: Dissimilarity Matrices**
#Spectral dis matrices
for (x in seq_along(extracted_name)){
  id <- extracted_name[[x]]
  name <- paste0("avg_pcs_scaled_", id)
  avg_pcs_scaled <- get(name)
  #significant change here. delete if it breaks something
  
  avg_pcs_scaled_sub <- avg_pcs_scaled %>% dplyr::select(PC1:PC4) 
  pcs_matrix <-stats::dist(avg_pcs_scaled_sub, method= "euclidean")
  pcs_matrix_name <- paste0("pcs_matrix_", id)
  assign(pcs_matrix_name, pcs_matrix)
}

#focal dis matrix
focal_matrix <- stats::dist(avg_focal_scaled, method = "euclidean")
#topo dis matrix
topo_matrix <- stats::dist(avg_topo_scaled, method = "euclidean")
#snow dis matrix:
snow_matrix <- stats::dist(avg_snow_scaled, method = "euclidean")








###**10 sets of mantels for each variable to determine ideal pcs**

Principal_Mantel <- function(a){ #Mantel of all spec versus other vars. "a" expects dissimilarity matrix to be compared with 1-10 pc spec dis matrices
all_results_df <- NULL
for(i in seq_along(extracted_name)) {
  id <- extracted_name[[i]]
  results_df <- data.frame( setNames(list(numeric(0)), paste0("R_", id)), setNames(list(numeric(0)), paste0("p_", id)))
  
  for(x in 1:10){
    name <- paste0("avg_pcs_scaled_", id)
    avg_pcs_scaled <- get(name)
    
    if(x == 1) {spec_sub <- avg_pcs_scaled[1]} else {spec_sub <- avg_pcs_scaled[,1:x]}
    spec_matrix <-stats::dist(spec_sub, method= "euclidean")
    
    
    spectral_mantel <- mantel(a, spec_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
    r <- spectral_mantel$statistic
    p <- spectral_mantel$signif
    row <- setNames(list(r, p), colnames(results_df))
    results_df<- rbind(results_df, row)
    
  }
  
  if(is.null(all_results_df)) {all_results_df <- results_df} else {all_results_df <- cbind(all_results_df, results_df)}
  
}
return(all_results_df)
}

veg_spec_results <- Principal_Mantel(focal_matrix)
enviro_spec_results <- Principal_Mantel(topo_matrix) #snow_matrix, topo_matrix, focal_matrix
snow_spec_results <- Principal_Mantel(snow_matrix)

veg_spec_results
enviro_spec_results
snow_spec_results

veg_spec_max <- c((which.max(veg_spec_results$R_003)), (which.max(veg_spec_results$R_016)), (which.max(veg_spec_results$R_018))) #list of three r_max values (003, 016, 018)
enviro_spec_max <- c((which.max(enviro_spec_results$R_003)), (which.max(enviro_spec_results$R_016)), (which.max(enviro_spec_results$R_018))) #list of three r_max values (003, 016, 018)
snow_spec_max <- c((which.max(snow_spec_results$R_003)), (which.max(snow_spec_results$R_016)), (which.max(snow_spec_results$R_018))) #list of three r_max values (003, 016, 018)

#export results to xl doc
write_xlsx(
  list(
    Sheet1 = veg_spec_results,
    Sheet2 = enviro_spec_results,
    Sheet3 = snow_spec_results
  ),
  path = "./Outputs/10_pcMantel_veg_enviro_snow_results.xls"
)



##**New spec distance matrices based on **##

#Spectral dis matrices


###**part 7: individual Mantel tests**
#veg spec mantels with appropriate pcs
for (x in seq_along(extracted_name)){
  id <- extracted_name[[x]]
    #generates new dis matrices with the appropriate #pcs for each test
    dis_name <- paste0("avg_pcs_scaled_", id)
    avg_pcs_scaled <- get(dis_name)
    max_pc <- veg_spec_max[x]
    max_pc_cols <- paste0("PC", 1:max_pc)
    
    avg_pcs_scaled_sub <- avg_pcs_scaled %>% dplyr::select(all_of(max_pc_cols)) 
    pcs_matrix <-stats::dist(avg_pcs_scaled_sub, method= "euclidean")
    pcs_matrix_name <- paste0("pcs_matrix_", id)
    assign(pcs_matrix_name, pcs_matrix)
 
  
  name <- paste0("pcs_matrix_", id)
  pcs_matrix <- get(name)
  
  grass_spectral_mantel <- mantel(focal_matrix, pcs_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
  grass_spectral_name <- paste0("grass_spectral_mantel_", id)
  assign(grass_spectral_name, grass_spectral_mantel)

} 
#enviro spec mantels with appropriate pcs
for (x in seq_along(extracted_name)){
  id <- extracted_name[[x]]
  #generates new dis matrices with the appropriate #pcs for each test
  dis_name <- paste0("avg_pcs_scaled_", id)
  avg_pcs_scaled <- get(dis_name)
  max_pc <- enviro_spec_max[x]
  max_pc_cols <- paste0("PC", 1:max_pc)
  
  avg_pcs_scaled_sub <- avg_pcs_scaled %>% dplyr::select(all_of(max_pc_cols)) 
  pcs_matrix <-stats::dist(avg_pcs_scaled_sub, method= "euclidean")
  pcs_matrix_name <- paste0("pcs_matrix_", id)
  assign(pcs_matrix_name, pcs_matrix)
  
  
  name <- paste0("pcs_matrix_", id)
  pcs_matrix <- get(name)
  
  topo_spectral_mantel <- mantel(topo_matrix, pcs_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
  topo_spectral_name <- paste0("topo_spectral_mantel_", id)
  assign(topo_spectral_name, topo_spectral_mantel)
  
} 
#snow spec mantels with appropriate pcs
for (x in seq_along(extracted_name)){
  id <- extracted_name[[x]]
  #generates new dis matrices with the appropriate #pcs for each test
  dis_name <- paste0("avg_pcs_scaled_", id)
  avg_pcs_scaled <- get(dis_name)
  max_pc <- snow_spec_max[x]
  max_pc_cols <- paste0("PC", 1:max_pc)
  
  avg_pcs_scaled_sub <- avg_pcs_scaled %>% dplyr::select(all_of(max_pc_cols)) 
  pcs_matrix <-stats::dist(avg_pcs_scaled_sub, method= "euclidean")
  pcs_matrix_name <- paste0("pcs_matrix_", id)
  assign(pcs_matrix_name, pcs_matrix)
  
  
  name <- paste0("pcs_matrix_", id)
  pcs_matrix <- get(name)
  
  snow_spectral_mantel <- mantel(snow_matrix, pcs_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)
  snow_spectral_name <- paste0("snow_spectral_mantel_", id)
  assign(snow_spectral_name, snow_spectral_mantel)
  
} 


topo_grass_mantel <- mantel(topo_matrix, focal_matrix, method ="spearman", permutations = 9999, na.rm = TRUE)

snow_grass_mantel <- mantel(snow_matrix, focal_matrix, method ="spearman", permutations = 9999, na.rm = TRUE)

topo_snow_mantel <- mantel(snow_matrix, topo_matrix,method ="spearman", permutations = 9999, na.rm = TRUE)
#could this significance be due to improper use of scaling?
#Results of all the tests:

#veg and pca distances by image
print(grass_spectral_mantel_003)
print(grass_spectral_mantel_016)
print(grass_spectral_mantel_018)

#environmental and pca distances by image
print(topo_spectral_mantel_003)
print(topo_spectral_mantel_016)
print(topo_spectral_mantel_018)

#snow and spec distances
print(snow_spectral_mantel_003)
print(snow_spectral_mantel_016)
print(snow_spectral_mantel_018)

#environmental and vegetation distances
print(topo_grass_mantel)

#snow and vegetation distances
print(snow_grass_mantel)

#snow~topo
#print(topo_snow_mantel) #fatal flaw with this analysis Ignore

results_grass_spectral <- c(grass_spectral_mantel_003$statistic, grass_spectral_mantel_016$statistic, grass_spectral_mantel_018$statistic)
base::plot(results_grass_spectral)
results_topo_spectral <- c(topo_spectral_mantel_003$statistic, topo_spectral_mantel_016$statistic, topo_spectral_mantel_018$statistic)
base::plot(results_topo_spectral)




##**Dan's Code**
f.multiplot <- function(l.plot, s.cols) {
  require(grid)
  s.numplots = length(l.plot)
  layout <- matrix(
    seq(1, s.cols * ceiling(s.numplots / s.cols)),
    ncol = s.cols,
    nrow = ceiling(s.numplots / s.cols)
  )
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
  for (i in 1:s.numplots) {
    matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
    print(l.plot[[i]], vp = viewport(layout.pos.row = matchidx$row, layout.pos.col = matchidx$col))
  }
} 

f.png <- function(st.file, v.dimuse) png(st.file, width = v.dimuse[1], height = v.dimuse[2], units = "in", res = v.dimuse[3], pointsize = v.dimuse[4])

Plot_Mantel <- function(s.use) ggplot(data = l.results[[s.use]], aes(x = 1:10, y = R_003)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = 4, color = "grey", linewidth = 0.75, linetype = "longdash") +
  geom_line(linewidth = 0.75, color = ch.col[1]) +
  geom_line(aes(y = R_016), linewidth = 0.75, color = ch.col[2]) +
  geom_line(aes(y = R_018), linewidth = 0.75, color = ch.col[3]) +
  geom_point(shape = -5*(l.results[[s.use]]$p_003 <= 0.050) + 21, size = 4, color = ch.col[1], fill = "white") +
  geom_point(aes(y = R_016), shape = -5*(l.results[[s.use]]$p_016 <= 0.050) + 21, size = 4, color = ch.col[2], fill = "white") +
  geom_point(aes(y = R_018), shape = -5*(l.results[[s.use]]$p_018 <= 0.050) + 21, size = 4, color = ch.col[3], fill = "white") +
  geom_point(data = df.legend, aes(x = xleg, y = (yleg + (s.use - 1)*100), color = lableg), size = 4) +
  geom_text(data = df.legend, aes(x = (xleg + 0.2), y = (yleg + (s.use - 1)*100), color = lableg, label = lableg), size = 5, hjust = 0) +
  scale_x_continuous(limits = c(0.5, 10.5), breaks = 1:10, labels = 1:10, expand = c(0, 0)) +
  scale_y_continuous(limits = c(-0.2, 0.6), breaks = seq(-0.2, 0.6, by = 0.1)) +
  scale_color_manual(values = ch.col) +
  xlab("Number of Components") +
  ylab("Mantel R-squared") +
  ggtitle(paste0("Spectral PCS ~", ch.results[s.use])) +
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.line = element_line(color = "black"),
    axis.ticks = element_line(color = "black"),
    plot.title = element_text(size = 16, color = "black", hjust = 0.5),
    panel.grid = element_blank(),
    panel.background = element_rect(fill = "white"),
    legend.position = "none"
  )


l.results <- list(veg_spec_results, enviro_spec_results, snow_spec_results)
ch.results <- c("Vegetation", "Environmental Variables", "Snow")
ch.col <- c("#1b9e77", "#d95f02", "#7570b3")
df.legend <- data.frame(xleg = rep(7.5, 3), yleg = c(-0.07, -0.12, -0.17), lableg = c("04/15/23", "10/09/24", "10/15/24"))

f.png("MantelSweeps.png", c(20, 6, 320, 15)) 
f.multiplot(lapply(1:3, Plot_Mantel), 3)
dev.off()


p1 <- Plot_Mantel(1)
p2 <- Plot_Mantel(2)
p3 <- Plot_Mantel(3)

# gridExtra
grid.arrange(p1, p2, p3, ncol = 3)







##**7.5: 4pc loadings**##

#for (x in seq_along(extracted_name[x])){ #loop not working for some reason
  id <- extracted_name[3] # change layers
  load_name <- paste0("loadings_", id)
  loadings <- get(load_name)
  loadings_df <- as.data.frame(as.matrix(loadings)[, 1:4])
  loadings_df <- loadings_df %>% mutate(across(everything(), abs)) #makes values absolute
  
  df_long <- loadings_df %>%
    mutate(bandname = rownames(.)) %>%
    mutate(bandname = as.numeric(bandname)) %>%
    pivot_longer(
      cols = -bandname,
      names_to = "Component",
      values_to = "Loading"
    )
  
  p <- ggplot(df_long, aes(x = bandname, y = Loading, color = Component, group = Component)) +
    geom_line() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    labs(x = "Band", y = "Loading") +
    theme_minimal()
  filename <- paste0(id, "loadings.jpg")
  ggsave(filename = file.path(output_folder, filename), plot = p,
         width = 10, height = 6, dpi = 300)
#}

  
  
  
  
###**Part 8: Plotting time**
#goals:
#I want to be able to be able to plug in two matrices into a function.
#I want that function to then apply the names of those matrices into the lables and output
######## Mantel test ########
scaled01 <- function(x){(x-min(x))/(max(x)-min(x))}#is this scaling redundant? 

for (x in seq_along(extracted_name)){
id <- extracted_name[[x]]
pcs_matrix <- get(paste0("pcs_matrix_", id))
spec<-scaled01(pcs_matrix) 
name <- paste0("spec_", id)
assign(name, spec)
}
veg<-scaled01(focal_matrix)
topo <- scaled01(topo_matrix)

#a: your y axis data, b: your x axis data
mantelrds <- function(a,b,id){

  mt <- mantel.randtest(a,b,nrepet=999)
  #### make data frame 
  VD <- as.data.frame(as.matrix(b))
  abbrev <- paste0("S_",alldata_avg$site)
  row.names(VD) <- abbrev
  colnames(VD) <- abbrev
  
  SD <- as.data.frame(as.matrix(a))
  abbrev <- paste0("S_",alldata_avg$site)
  row.names(SD) <- abbrev
  colnames(SD) <- abbrev
  
  ##### combine in one dataframe 
  colnames(VD)==colnames(SD)
  row.names(VD)==row.names(SD)
  
  nami <- rep(colnames(VD),each=ncol(VD))
  VVD <- as.vector(as.matrix(VD))
  VSD <- as.vector(as.matrix(SD))
  
  dat <- as.data.frame(cbind(VVD, VSD))
  dat <- cbind(nami, dat)
  
  
  plotID <- rep(abbrev, times=ncol(VD))
  VD_SD_frame <- cbind(dat, plotID)
  
  rds_name_1 <- paste0("./Outputs/mantel_rds/", id, "_res_mantel_sites.rds")
  saveRDS(mt, file=rds_name_1)
  rds_name_2 <- paste0("./Outputs/mantel_rds/" ,id,"_VD_SD_frames.rds")
  saveRDS(VD_SD_frame, file=rds_name_2)
  
  ###axis labels
  x_axis <- deparse(substitute(b))
  if (deparse(substitute(a)) == "Environmental Distances") {
    y_axis <- "Environmental Distances"
  } else {
    y_axis <- spec_date
  }
  plot_title <-paste0(x_axis, "~", y_axis) 
  #######################
  #### Plots per Site
  
  VD_SD_frame <- readRDS(rds_name_2) 
  # coliix <- dichromat(colours = rainbow(n=ncol(MFD)),type="deutan")
  # coliix <- dichromat(colours = rainbow(n=ncol(MFD)),type="protan")
  # plot(1:length(coliix), 1:length(coliix), col=coliix, pch=19, cex=3, xlab="", ylab="")
  
  
  dat <- dat %>% mutate(siteID=substr(nami,1,4)) %>% group_by(siteID)
  
  coliix <- dichromat(colours = rainbow(n=max(table(dat$nami))),type="tritan")
  

  #### one Site
  ggplot(data=dat, aes(x=VVD, y=VSD)) +
    scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
    scale_x_continuous(breaks=c(0,0.5,1))+
    geom_smooth(aes(color=plotID),method="lm",lwd=0.6,
                show.legend = T, se = F, fullrange=F)+
    geom_point(color="grey30", pch=1, size=2)+
    # scale_linetype_manual(values=c(lini))+
    scale_color_manual(values=coliix)+
    # scale_color_viridis_d(option="viridis")+
    labs(x=x_axis,
         y=y_axis,
         title = plot_title)+
    # geom_smooth(method = "lm", formula = y~-1+x+I(x^2), col="black", lwd=0.7,
    #             lty=2, alpha=0.4, fullrange=T)
    geom_smooth(method = "lm", color="black",
                formula = y~x,  lwd=1.2,
                lty=1, alpha=0.4, fullrange=T, level=0.95)+
    # ggtitle(label = substr(namFD[1],1,4))+
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey",size=0.2),
          # strip.text= element_text(size=12),
          # legend.text = element_text(size = 9),
          # legend.title = element_text(size = 12),
          # legend.box.spacing = unit(1.8,"line"),
          # legend.key.width = unit(2.2,"line"),
          axis.title.x = element_text(margin = margin(t = 15), size=18),
          axis.title.y = element_text(margin = margin(r = 15), size=18),
          plot.title =element_text(size=14, face="bold", hjust = -0.15),
          legend.key = element_rect(size = 0.1, fill = "white"),
          legend.key.size = unit(1.2,"line"),
          # text=element_text(size=12),
          axis.ticks = element_line(colour = "grey", size = 0.1)) 
  ggsaveid_1 = paste0("./Outputs/single_site",id,".jpg")
  ggsave(ggsaveid_1) 
  
  
  
  #################
  ##### STRIP PLOT
  ggplot(data=dat, aes(x=VVD, y=VSD)) +
    scale_y_continuous(limit= c(0,1),breaks=c(0,0.5,1))+
    scale_x_continuous(breaks=c(0,0.5,1))+
    geom_point(color="grey60", pch=1, size=1.2)+
    geom_smooth(aes(group = plotID),method="lm",color="grey10" ,lwd=0.4,
                show.legend = F, se = F, fullrange=F)+
    facet_wrap(~siteID,scales = "free")+
    # scale_linetype_manual(values=c(lini))
    # scale_color_manual(values=coliix)+
    # scale_color_viridis_d(option="viridis")+
    labs(x=x_axis,y=y_axis)+
    geom_smooth(method = "lm", formula = y~x, col="red3", lwd=1,
                lty=1, alpha=0.4, fullrange=T, level=0.95)+
    # ggtitle(label = substr(namFD[1],1,4))+
    theme(panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.grid.major =  element_line(colour = "grey",size=0.2),
          strip.text= element_text(size=14),
          legend.text = element_text(size = 9),
          legend.title = element_text(size = 12),
          legend.box.spacing = unit(1.8,"line"),
          legend.key.width = unit(2.2,"line"),
          axis.title.x = element_text(margin = margin(t = 15), size=18),
          axis.title.y = element_text(margin = margin(r = 15), size=18),
          plot.title =element_text(size=14, face="bold", hjust = -0.15),
          legend.key = element_rect(size = 0.1, fill = "white"),
          legend.key.size = unit(1.2,"line"),
          text=element_text(size=12),
          axis.ticks = element_line(colour = "grey", size = 0.1)) 
  ggsaveid_2 <- paste0("./Outputs/strip_plot",id,".jpg")
  ggsave(ggsaveid_2,width = 12, height = 8)

}

#run function 
for( x in seq_along(extracted_name)){
  if(extracted_name[[x]] == "003") {
    spec_date <- "Spectral Componant Distances 003"} else {if(extracted_name[[x]] == "016") {spec_date <- "Spectral Componant Distances 016"} else {spec_date <- "Spectral Componant Distances 018"}}
spec_id <- paste0("spec_",extracted_name[[x]])
name <- get(spec_id)

mantelrds(name, veg, id = paste0(spec_date, "~ Focal Grass Species Abundance Distances"))
mantelrds(name, topo, id = paste0(spec_date, "~ Environmental Distances") )
}

mantelrds(veg, topo, id = "Focal Grass Species Abundance Distances ~ Environmental Distances")




