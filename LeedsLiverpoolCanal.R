# ============================================================
# Leeds-Liverpool Canal LiDAR DSM + RGB tiles -> 3D GLB
# ROBUST LARGE-MOSAIC VERSION
# ============================================================

# Packages
library(terra)
library(rayshader)
library(rgl)
library(rgl2gltf)

# ------------------------------------------------------------
# USER INPUTS
# ------------------------------------------------------------

dsm_file <- "LeedsLiverpoolCanal_LidarDSM/Download_2979462/england-dsm-fr-1m_6392693/se/se23se_fz_dsm_1m.tif"

rgb_files <- c(
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2633_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2634_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2733_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2734_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2833_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2834_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2933_rgb_250_03.jpg",
  "LeedsLiverpoolCanal_Aerial/Download_2979463/getmapping_rgb_25cm_6392694/se/se2934_rgb_250_03.jpg"
)

out_glb <- "LeedsLiverpoolCanal_dsm.glb"

# ------------------------------------------------------------
# SETTINGS
# ------------------------------------------------------------

# Increase for smaller/more stable mesh
geom_res_m <- 2.5

# Reduce texture size if needed
max_texture_dim <- 8000

# Vertical exaggeration
zscale_val <- 2

# Camera
theta_val <- 135
phi_val   <- 35
zoom_val  <- 0.9
fov_val   <- 0

# Texture shading
shading_mode <- "rgb_only"

# ------------------------------------------------------------
# 1) LOAD DSM
# ------------------------------------------------------------

dsm <- rast(dsm_file)

# Fill NA areas
dsm[is.na(dsm)] <- 3

# ------------------------------------------------------------
# 2) RGB HELPERS
# ------------------------------------------------------------

to_rgb3 <- function(r){
  
  nl <- nlyr(r)
  
  if(nl == 4){
    r <- r[[1:3]]
  }
  
  if(nl == 1){
    r <- c(r,r,r)
  }
  
  if(nlyr(r) != 3){
    stop("RGB raster does not have 3 bands")
  }
  
  names(r) <- c("R","G","B")
  
  r
}

read_rgb <- function(f, target_crs){
  
  r <- rast(f)
  
  # Assign DSM CRS if missing
  if(is.na(crs(r)) || crs(r) == ""){
    crs(r) <- target_crs
  }
  
  to_rgb3(r)
}

# ------------------------------------------------------------
# 3) LOAD + MOSAIC RGB
# ------------------------------------------------------------

rgb_list <- lapply(
  rgb_files,
  read_rgb,
  target_crs = crs(dsm)
)

rgb_mosaic <- do.call(mosaic, rgb_list)

# ------------------------------------------------------------
# 4) FIND OVERLAP
# ------------------------------------------------------------

overlap_ext <- intersect(
  ext(dsm),
  ext(rgb_mosaic)
)

dsm_crop <- crop(
  dsm,
  overlap_ext,
  snap = "in"
)

rgb_crop <- crop(
  rgb_mosaic,
  overlap_ext,
  snap = "in"
)

# ------------------------------------------------------------
# 5) RESAMPLE RGB TO DSM GRID
# ------------------------------------------------------------

rgb_resampled <- resample(
  rgb_crop,
  dsm_crop,
  method = "bilinear"
)

# ------------------------------------------------------------
# 6) MASK VALID AREAS
# ------------------------------------------------------------

valid_mask <- app(
  rgb_resampled,
  fun = function(x) any(!is.na(x))
)

dsm_masked <- mask(
  dsm_crop,
  valid_mask,
  updatevalue = NA
)

rgb_masked <- mask(
  rgb_resampled,
  valid_mask,
  updatevalue = NA
)

dsm_masked[is.na(dsm_masked)] <- 3

# ------------------------------------------------------------
# 7) TRIM
# ------------------------------------------------------------

dsm_masked <- trim(dsm_masked)
rgb_masked <- trim(rgb_masked)

common_ext <- intersect(
  ext(dsm_masked),
  ext(rgb_masked)
)

dsm_masked <- crop(dsm_masked, common_ext)
rgb_masked <- crop(rgb_masked, common_ext)

# ------------------------------------------------------------
# 8) OPTIONAL GEOMETRY DOWNSAMPLING
# ------------------------------------------------------------

geom_template <- rast(
  ext = ext(dsm_masked),
  crs = crs(dsm_masked),
  resolution = geom_res_m
)

dsm_geom <- resample(
  dsm_masked,
  geom_template,
  method = "bilinear"
)

rgb_geom <- resample(
  rgb_masked,
  geom_template,
  method = "bilinear"
)

# ------------------------------------------------------------
# 9) FORCE IDENTICAL GRID
# ------------------------------------------------------------

rgb_geom <- crop(
  rgb_geom,
  dsm_geom,
  snap = "in"
)

rgb_geom <- extend(
  rgb_geom,
  dsm_geom
)

rgb_geom <- resample(
  rgb_geom,
  dsm_geom,
  method = "bilinear"
)

# ------------------------------------------------------------
# 10) LIMIT MASSIVE TEXTURES
# ------------------------------------------------------------

if(
  ncol(rgb_geom) > max_texture_dim ||
  nrow(rgb_geom) > max_texture_dim
){
  
  scale_factor <- max(
    ncol(rgb_geom) / max_texture_dim,
    nrow(rgb_geom) / max_texture_dim
  )
  
  new_res <- res(rgb_geom) * scale_factor
  
  tex_template <- rast(
    ext = ext(rgb_geom),
    crs = crs(rgb_geom),
    res = new_res
  )
  
  rgb_geom <- resample(
    rgb_geom,
    tex_template,
    method = "bilinear"
  )
  
  dsm_geom <- resample(
    dsm_geom,
    tex_template,
    method = "bilinear"
  )
}

# ------------------------------------------------------------
# 11) SAVE QA PNG
# ------------------------------------------------------------

png(
  "LeedsLiverpoolCanal_rgb_fullres.png",
  width = ncol(rgb_geom),
  height = nrow(rgb_geom)
)

plotRGB(rgb_geom, scale = 255)

contour(
  dsm_geom,
  add = TRUE,
  col = "white",
  lwd = 0.5
)

dev.off()

# ------------------------------------------------------------
# 12) BUILD RAYSHADER MATRICES
# ------------------------------------------------------------

elmat <- raster_to_matrix(dsm_geom)

nr <- nrow(dsm_geom)
nc <- ncol(dsm_geom)

rgb_arr <- array(
  NA_real_,
  dim = c(nr, nc, 3)
)

for(i in 1:3){
  
  vals <- values(
    rgb_geom[[i]],
    mat = FALSE
  )
  
  rgb_arr[,,i] <- matrix(
    vals,
    nrow = nr,
    ncol = nc,
    byrow = TRUE
  ) / 255
}

# Flip vertically for rayshader
# rgb_arr <- rgb_arr[nr:1,,]

# ------------------------------------------------------------
# 13) BUILD TEXTURE
# ------------------------------------------------------------

if(shading_mode == "rgb_only"){
  
  tex <- rgb_arr
  
} else if(shading_mode == "ambient"){
  
  amb <- ambient_shade(
    elmat,
    zscale = zscale_val
  )
  
  tex <- add_overlay(
    rgb_arr,
    amb,
    alphalayer = 0.12
  )
  
} else {
  
  ray <- ray_shade(
    elmat,
    sunaltitude = 60
  )
  
  tex <- add_overlay(
    rgb_arr,
    ray,
    alphalayer = 0.10
  )
}

# ------------------------------------------------------------
# 14) PREVIEW
# ------------------------------------------------------------

rgl::clear3d()

plot_3d(
  hillshade = tex,
  heightmap = elmat,
  zscale    = zscale_val,
  theta     = theta_val,
  phi       = phi_val,
  zoom      = zoom_val,
  fov       = fov_val
)

# ------------------------------------------------------------
# 15) SAVE PREVIEW
# ------------------------------------------------------------

png(
  "LeedsLiverpoolCanal_dsm_preview.png",
  width = 2000,
  height = 2000
)

plot_map(tex)

dev.off()

# ------------------------------------------------------------
# 16) EXPORT GLB
# ------------------------------------------------------------

scene <- scene3d()

gltf <- rgl2gltf::as.gltf(scene)

rgl2gltf::writeGLB(
  gltf,
  out_glb
)

cat("\nGLB export complete:\n")
cat(out_glb, "\n")