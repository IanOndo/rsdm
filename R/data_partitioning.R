#'                                 Geographic "block" partitioning of occurrence records for evaluation
#'
#' Partition both occurrence records and background into evaluation bins based on some spatial rules.
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param bg A RasterLayer object of the background to be partitionned.
#' @param assign_method A character string specifying the method to use to assign background cell to each block between: `nearest.cluster` and `nearest.center` (the default).
#' @param grid_res A numeric value specifying the resolution of the output raster layer in decimal degree. Default is 0.1666667 (~ 20 km).
#' @param coastline An optional character string specifying the path to a shapefile delineating the coastline.
#' @param sf A logical. Should a sf object be returned ? Default is TRUE
#' @param ..., Additional parameters to be passed to \code{pamEqual or kmeansEqual} or \code{make_geographic_domain}
#' @return A sf object (default) or a RasterStack of the geographically masked blocks of occurrence and background points/cells.
#' @export
make_geographic_block <- function(kdat, k = 3, bg = NULL, assign_method=c("nearest.center"), grid_res = 0.1666667, coastline = NULL, sf=TRUE, verbose=TRUE, ...){

  #==================================
  # 1. Partition occurrence records #
  #==================================
  if(verbose) cat('> [...partitioning occurrence records in ',k,' groups...]\n')

  list_args_partalgo<- list(kdat=kdat, k=k, verbose=FALSE, plot=FALSE)
  dot.args <- list(...)
  if(length(dot.args)>0L){
    partalgo.args <- dot.args[names(dot.args) %in% c("iter_max", "method", "nstart", "n.outliers")]
    domain.args <- dot.args[!names(dot.args)%in%names(partalgo.args)]
    if(length(partalgo.args)>0L){
      list_args_partalgo<- append(list_args_partalgo, partalgo.args)
    }
  }

  kdatAssigned <- try(do.call(pamEqual,list_args_partalgo),silent=TRUE)

  # if it fails, try with a different algorithm
  if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
    alt_method <- if(!is.null(list_args_partalgo$method)) setdiff(c("utmzone.outliers","top.outliers","iqr.outliers"),list_args_partalgo$method) else "utmzone.outliers"
    list_args_partalgo$method <- alt_method
    kdatAssigned <- try(do.call(pamEqual,list_args_partalgo),silent=TRUE)

    if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
      alt_method <- setdiff(c("utmzone.outliers","top.outliers"), alt_method)
      list_args_partalgo$method <- alt_method
      kdatAssigned <- try(do.call(pamEqual,list_args_partalgo),silent=TRUE)

      if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
        if(inherits(kdatAssigned,"try-error")) print(kdatAssigned) else cat("PAM partioning algorithm did not converge\n")

        if(verbose) cat("Trying kmeans algorithm...\n")

        kdatAssigned <- try(do.call(kmeansEqual,list_args_partalgo),silent=TRUE)

        # if it fails, try with a different algorithm
        if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
          alt_method <- if(!is.null(list_args_partalgo$method)) setdiff(c("utmzone.outliers","top.outliers","iqr.outliers"),list_args_partalgo$method) else "utmzone.outliers"
          list_args_partalgo$method <- alt_method
          kdatAssigned <- try(do.call(kmeansEqual,list_args_partalgo),silent=TRUE)

          if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
            alt_method <- setdiff(c("utmzone.outliers","top.outliers"), alt_method)
            list_args_partalgo$method <- alt_method
            kdatAssigned <- try(do.call(kmeansEqual,list_args_partalgo),silent=TRUE)

            if(inherits(kdatAssigned,"try-error") || !kdatAssigned$converged){
              if(inherits(kdatAssigned,"try-error")) print(kdatAssigned) else cat("Kmeans partioning algorithm did not converge\n")
              stop("Unable to partition the data")
            }
          }

        }

      }
    }

  }

  # make sure each group has at least 2 points
  kdatAssign <- kdatAssigned
  process=any(table(kdatAssign$Data$assigned)==1) && list_args_partalgo$k>2
  while(process){

    if(verbose) cat(">>...making sure to keep >=2 points in each geographic blocks...\n")

    list_args_partalgo$k <- list_args_partalgo$k-1

    kdatAssign <- try(do.call(pamEqual,list_args_partalgo),silent=TRUE)

    if(inherits(kdatAssign,"try-error"))
      kdatAssign <- try(do.call(kmeansEqual,list_args_partalgo),silent=TRUE)

    process = !(inherits(kdatAssign,"try-error") || !kdatAssign$converged) && any(table(kdatAssign$Data$assigned)==1) && list_args_partalgo$k>2
  }

  if(inherits(kdatAssign,"try-error")){
    print(kdatAssign)
    stop("Unable to partition the data")
  }else if(!kdatAssign$converged){
    cat("Partioning algorithm did not converge with minimum 2 points per geographic block\n")
  }else{
    kdatAssigned <- kdatAssign
    if(list_args_partalgo$k!=k){
      if(verbose) cat(">>...changing the number of geographic block to:",list_args_partalgo$k,"\n")
      k <- list_args_partalgo$k
    }
  }

  if(any(table(kdatAssigned$Data$assigned)==1))
    warning(">[Partition of occurrence records]: at least one group contains 1 unique record...<")

  #==========================================
  # 2. Build convex hulls around each group #
  #==========================================
  if(verbose) cat('> [...building convex hulls around groups...]\n')
  proj <- if(is.null(bg) || is.na(sf::st_crs(bg))) '+proj=longlat +datum=WGS84' else sf::st_crs(bg)
  sf_loc_data <- kdatAssigned$Data %>%
    sf::st_as_sf(coords=c('lon', 'lat'), crs=sf::st_crs(proj))
  # create the first convex hull
  j = 1
  sf_loc_grp <- sf_loc_data %>%
    dplyr::filter(assigned==j)

  sf_block <- sf::st_as_sf(data.frame(kgroup=j, geometry=sf::st_buffer(sf::st_convex_hull(sf::st_union(sf_loc_grp)), dist=0.01)))

  # add other hulls
  j = j + 1
  while(j <= k){
    sf_loc_grp <- sf_loc_data %>%
      dplyr::filter(assigned==j)

    sf_block <- sf_block %>%
      rbind(sf::st_sf(data.frame(kgroup=j, geometry=sf::st_buffer(sf::st_convex_hull(sf::st_union(sf_loc_grp)), dist=0.01))))
    j = j + 1
  }
  #=============================
  # 3. Build disjoint polygons #
  #=============================
  sf_block_disjoint <- sf_block
  # compute intersections between all polygons and retrieve polygons of the intersections
  sf_block_intersect <- sf::st_intersection(sf_block) %>%
    dplyr::filter(n.overlaps>1) %>%
    dplyr::select(kgroup)
  # substract polygons intersection to the blocks to obtain block disjoints
  if(nrow(sf_block_intersect)>0L){
    sf_block_disjoint <- suppressMessages(sf_block_intersect %>%
      sf::st_difference(sf_block, .) %>%
      dplyr::select(kgroup))
  }
  # if(!all(sf::st_is(sf_block_disjoint ,"POLYGON"))){
  #   list_simple_features <- 1:length(sf::st_geometry(sf_block_disjoint)) %>% lapply(function(p) sf_block_disjoint[p,] %>% sf::st_cast("POLYGON",warn=FALSE))
  #   sf_block_disjoint <- do.call(rbind, list_simple_features)
  # }

  #=========================
  # 4. Make the background #
  #=========================
  if(verbose) cat('> [...making the background...]\n')
  if(is.null(bg)){

    if(is.null(coastline))
      coastline = RGeodata::terrestrial_lands

    domain <- try(do.call(make_geographic_domain, append(list(loc_dat=kdat, verbose=FALSE, land_file=coastline), domain.args)), silent=TRUE)

    if(inherits(domain,'try-error')){
      cat(domain)
      stop("Unable to build a background from the set of points")
    }

    # fasterize the domain
    raster_template  <- tryCatch(raster::raster(domain, res = grid_res), error=function(err) raster::raster(as(domain,"Spatial"), res = grid_res))
    bg <- suppressMessages(fasterize::fasterize(sf::st_cast(domain,"MULTIPOLYGON"), raster_template) - 1 )# -1 to ensure background is a distinct group
  }else if(!inherits(bg,"RasterLayer")){
    stop('bg must be a RasterLayer object')
  }
  #----------------------------------------------------
  #= 4.a Assign background cells to the block cluster #
  #----------------------------------------------------
  if(verbose) cat('     > [...assign background cells to the block cluster...]\n')
  #bg_fasterized <- raster::mask(fasterize::fasterize(sf_block_disjoint, bg, field="kgroup", background=0), bg)
  # make sure all groups are composed of polygons
  sf_block_disjoint %<>% sf::st_collection_extract(type="POLYGON" )%>% sf::st_make_valid()
  bg_fasterized <- raster::mask(fasterize::fasterize(sf::st_cast(sf_block_disjoint,"POLYGON",warn=FALSE), bg, field="kgroup", background=0), bg)
  # assign occurrence records cells outside of the block to their corresponding group
  # if(nrow(sf_block_intersect)>0L){
  #   cell_occ <- kdatAssigned$Data %>%
  #     sf::st_as_sf(coords=c('lon', 'lat'), crs=sf::st_crs(bg)) %>%
  #     sf::st_intersects(., sf_block_disjoint) %>%
  #     purrr::map_lgl(., function(x) length(x) == 0L ) %>%
  #     dplyr::filter(kdatAssigned$Data, .)
  # }
  tbl <- data.table::data.table(cell=which(raster::values(bg_fasterized)>=0), key="cell")
  tbl[,`:=`(Group = as.character(bg_fasterized[cell]),
            Lon = raster::xFromCell(bg_fasterized, cell),
            Lat = raster::yFromCell(bg_fasterized, cell))]
  data.table::setkey(tbl, 'Group')

  switch(assign_method,
       "nearest.center" = {
         centers <- kdatAssigned$Centers
         tbl[, kgroup := ifelse(Group=='0', as.character(which.min(geosphere::distGeo(cbind(Lon, Lat), centers))), Group), by=seq_len(NROW(tbl))]
         bg_fasterized[tbl['0'][['cell']]] <- tbl['0'][['kgroup']]
       },
       "nearest.cluster" = {
         kgroup <- tbl['0'][, .(Lon,Lat)] %>%
           sf::st_as_sf(coords=c('Lon', 'Lat'), crs=sf::st_crs(bg)) %>%
           sf::st_distance(., sf_block_disjoint) %>%
           min.col() %>% as.character()
         tbl[, kgroup := replace(Group, Group=='0', kgroup)]
         bg_fasterized[tbl['0'][['cell']]] <- as.numeric(tbl['0'][['kgroup']])
       }
  )

  #-------------------------------
  #= 4.b Create background masks #
  #-------------------------------
  if(verbose) cat('     > [...create background masks...]\n')
  data.table::setkey(tbl, 'kgroup')
  comb <- utils::combn(1:k,k-1)
  if(!sf){
    bg_masks <- stack()
    for(j in 1:k){
      blcks <- comb[,j]
      blck <- setdiff(1:k,blcks)
      cells_group <- subset(tbl, kgroup %in% blcks, select=cell)$cell#c(tbl[as.character(blcks[1])][['cell']], tbl[as.character(blcks[2])][['cell']])
      non_cells_group <- tbl[as.character(blck)][['cell']]
      bg_group <- bg_fasterized
      bg_group[cells_group] <- "1"
      bg_group[non_cells_group] <- "0"
      bg_masks <- raster::addLayer(bg_masks, bg_group)
    }
    names(bg_masks) <- apply(comb,2,function(x) paste(letters[x],collapse="_"))
  }else{
    bg_masks <- bg_fasterized %>%
      raster::setValues(as.numeric(raster::values(.))) %>%
      stars::st_as_stars() %>%
      sf::st_as_sf(as_points=FALSE, merge=TRUE) %>%
      sf::st_cast("MULTIPOLYGON") %>%
      sf::st_make_valid() %>%
      st_remove_holes()

    # switch off Spherical geometry (s2)
    if(any(!bg_masks %>% sf::st_is_valid())){
      sf_use_s2(FALSE)
      bg_masks <- bg_fasterized %>%
        raster::setValues(as.numeric(raster::values(.))) %>%
        stars::st_as_stars() %>%
        sf::st_as_sf(as_points=FALSE, merge=TRUE) %>%
        sf::st_cast("MULTIPOLYGON") %>%
        sf::st_make_valid()
      sf_use_s2(TRUE)
    }
  }
  if(verbose) cat('> [...completed...]\n')
  # return multi-layer background mask
  return(bg_masks)
}

#'                                 Kmeans clustering with near equal group size
#'
#' Modified from https://rviews.rstudio.com/2019/06/13/equal-size-kmeans/
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param iter_max A numeric integer specifying the maximum number of iterations to run.
#' @return A list with two elements:
#'         \code{Data} A two-column matrix or data.frame or data.table with longitude and latitude coordinates and a column 'assigned' specifying
#'         which of the k cluster each point belongs to.
#'         \code{Centers} A two-column data.frame with longitude and latitude coordinates of the cluster centers.
#' @export
kmeansEqual <- function(kdat, k = 3, iter_max = 30, tolerance.iter=1e-08, random.seed=1234, method=c("iqr.outliers","top.outliers","utmzone.outliers"), nstart=20, n.outliers=10, verbose=TRUE, plot=TRUE){

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 0: Check input data\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  if(!any(inherits(kdat, c("matrix","data.frame","data.table"))))
    stop("Argument kdat must be a matrix, a data.frame or a data.table")
  if(nrow(kdat) < 2)
    stop("Input data must have at least 2 coordinates.")

  if(!is.numeric(k))
    stop("Argument 'k' must be numeric.")
  if(k<0)
    stop("The number of cluster must be > 0.")
  if(k<2)
    warning("The number of group requested is < 2.")

  if(!is.numeric(iter_max))
    stop("Argument 'iter_max' must be numeric")
  if(iter_max<0)
    stop("The number of iterations must be > 0.")

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 1: Initialize "naive" clusters with a basic kmeans\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  set.seed(random.seed)

  kclust <- kdat %>%
    stats::kmeans(k, nstart = nstart)

  # initial clusters
  centers = kclust$centers

  # transform to cartesian coordinates
  centers %<>%
    as.data.frame(row.names=NULL) %>%
    sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
    sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
    sf::st_coordinates()

  if(is.matrix(kdat))
    kdat %<>%
    as.data.frame(row.names=NULL)

  # retrieve coordinates if needed
  if(!all(c("lon","lat") %in% colnames(kdat))){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(kdat))[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(kdat))[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
    kdat %<>%
      dplyr::rename(lon = coordHeaders[1], lat = coordHeaders[2])
  }

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 2: Re-assign point to clusters\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  # start iterations
  iter = 0
  process = TRUE
  while(iter <= iter_max & process){

    if(verbose) cat('iteration: ',iter+1,'\n')

    kclust <- kdat %>%
      dplyr::select(lon,lat) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
      sf::st_coordinates() %>%
      stats::kmeans(centers)

    kdat$assigned = kclust$cluster
    #---------------------------------------------------------------------------
    #=step 2: Compute the distance matrix between each point and the centroids
    #---------------------------------------------------------------------------

    # compute mean distance between centers
    distMean <- mean(stats::dist(centers))

    # transform back to lon/lat
    centers %<>%
      as.data.frame(row.names=NULL) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
      sf::st_transform(crs=sf::st_crs(4326)) %>%
      sf::st_coordinates()

    # calculate the distance to the center
    for(j in 1:nrow(centers)){
      kdat <- kdat %<>%
        dplyr::mutate(kdist = geosphere::distGeo(cbind(lon,lat), centers[j,]))
      # rename the column
      kdat <- kdat %>%
        dplyr::rename(!!paste0('kdist',j) := kdist)
    }
    if(iter > 0){
      # assign the outliers from each group to the closest cluster center
      switch(method[1],

             "top.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                   )
               }
             },

             "iqr.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                  .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                            .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             },

             "utmzone.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             }
      )
      kdat.copy %<>% data.frame(row.names=NULL)
      kdat <- kdat.copy

    }
    #---------------------------------------------------------------------------
    #=step 3: Re-assign points to clusters
    #---------------------------------------------------------------------------
    kdat$index = 1:nrow(kdat)
    working = kdat
    ReAssignment = nrow(kdat) - (nrow(kdat) %% k)

    for(i in 1:ReAssignment){
      #cluster counts can be off by 1 due to uneven multiples of k.
      j = if(i %% k == 0) k else (i %% k)
      itemloc =  working$index[which.min(working[,(paste0("kdist", j))])][1]
      kdat$assigned[kdat$index == itemloc] = j
      working %<>%
        dplyr::filter(!index == itemloc)
    }
    # if there are some leftover points, assign to whoever's closest, without regard to k
    if(length(working$index)>0L){
      for(i in working$index){
        kdat$assigned[kdat$index == i] = which.min(working[working$index == i,grepl("kdist",colnames(working))])
      }
    }

    # assign the outliers from each group to the closest cluster center
    switch(method[1],

           "top.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                 )
             }
           },

           "iqr.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                          .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           },

           "utmzone.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           }
    )
    kdat.copy %<>% data.frame(row.names=NULL)
    kdat <- kdat.copy
    if(verbose){
      cat('#---------------------------------------------------------------------------\n')
      cat('#=step 4: Recalculate the centroids\n')
      cat('#---------------------------------------------------------------------------\n')
    }
    j = 1
    NewCenters <- kdat %>%
      dplyr::filter(assigned == j) %>%
      dplyr::select(lon, lat) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
      sf::st_coordinates() %>%
      stats::kmeans(1) %$% centers
    while(j < k){
      j = j+1
      NewCenters %<>%
        rbind(kdat %>%
                dplyr::filter(assigned == j) %>%
                dplyr::select(lon, lat) %>%
                sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
                sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
                sf::st_coordinates() %>%
                stats::kmeans(1) %$% centers)
    }
    # New centroids
    NewCenters %<>% data.frame(row.names=NULL)
    centers <- NewCenters

    # keep assigments lon, lat and assigments only
    kdat <- kdat %>%
      dplyr::select(lon, lat, assigned)

    if(plot){
      cent <- centers %>%
        as.data.frame(row.names=NULL) %>%
        sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
        sf::st_transform(crs=sf::st_crs(4326)) %>%
        sf::st_coordinates() %>%
        as.data.frame(row.names=NULL)

      if(!is.factor(kdat$assigned))
        kdat$assigned %<>% as.factor()
      print(
        kdat %>% ggplot2::ggplot(ggplot2::aes(x = lon, y = lat, color = assigned)) +
          ggplot2::theme_minimal() + ggplot2::geom_point()  +
          ggplot2::geom_point(data = cent, ggplot2::aes(x = X, y = Y), color = "black", size = 4)
      )
    }

    process = abs(mean(stats::dist(centers)) - distMean) > tolerance.iter
    iter = iter + 1
  }
  # end iterations

  # assignment to factor
  if(!is.factor(kdat$assigned))
    kdat$assigned %<>% as.factor()

  centers %<>%
    sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
    sf::st_transform(crs=sf::st_crs(4326)) %>%
    sf::st_coordinates()

  return(list(Data=kdat, Centers=centers, converged=!process))
}

#'                                 Pam clustering with near equal group size
#'
#'
#'
#' @param kdat A two-column matrix or data.frame or data.table with longitude and latitude coordinates.
#' @param k A numeric integer specifying the number of cluster required.
#' @param iter_max A numeric integer specifying the maximum number of iterations to run.
#' @return A list with two elements:
#'         \code{Data} A two-column matrix or data.frame or data.table with longitude and latitude coordinates and a column 'assigned' specifying
#'         which of the k cluster each point belongs to.
#'         \code{Centers} A two-column data.frame with longitude and latitude coordinates of the cluster centers.
#' @export
pamEqual <- function(kdat, k = 3, iter_max = 30, tolerance.iter=1e-08, random.seed=1234, method=c("iqr.outliers","top.outliers","utmzone.outliers"), nstart=20, n.outliers=10, verbose=TRUE, plot=TRUE){

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 0: Check input data\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  if(!any(inherits(kdat, c("matrix","data.frame","data.table"))))
    stop("Argument kdat must be a matrix, a data.frame or a data.table")
  if(nrow(kdat) < 2)
    stop("Input data must have at least 2 coordinates.")

  if(!is.numeric(k))
    stop("Argument 'k' must be numeric.")
  if(k<0)
    stop("The number of cluster must be > 0.")
  if(k<2)
    warning("The number of group requested is < 2.")

  if(!is.numeric(iter_max))
    stop("Argument 'iter_max' must be numeric")
  if(iter_max<0)
    stop("The number of iterations must be > 0.")

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 1: Initialize "naive" clusters with a basic pam\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  set.seed(random.seed)

  kclust <- kdat %>%
    cluster::pam(k, nstart = nstart)


  # initial clusters
  centers = kclust$medoids

  # transform to cartesian coordinates
  centers %<>%
    as.data.frame(row.names=NULL) %>%
    sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
    sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
    sf::st_coordinates()

  if(is.matrix(kdat))
    kdat %<>%
    as.data.frame(row.names=NULL)

  # retrieve coordinates if needed
  if(!all(c("lon","lat") %in% colnames(kdat))){
    id_x_lon 	<- grep(pattern = "[Ll][Oo][Nn]|[Xx]",x = names(kdat))[1]
    id_y_lat 	<- grep(pattern = "[Ll][Aa][Tt]|[Yy]",x = names(kdat))[1]
    coordHeaders <- c(id_x_lon, id_y_lat)
    kdat %<>%
      dplyr::rename(lon = coordHeaders[1], lat = coordHeaders[2])
  }

  if(verbose){
    cat('#---------------------------------------------------------------------------\n')
    cat('#=step 2: Re-assign point to clusters\n')
    cat('#---------------------------------------------------------------------------\n')
  }

  # start iterations
  id.centers = kclust$id.med
  iter = 0
  process = TRUE
  while(iter <= iter_max & process){

    if(verbose) cat('iteration: ',iter+1,'\n')

    kclust <- kdat %>%
      dplyr::select(lon,lat) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
      sf::st_coordinates() %>%
      cluster::pam(k,medoids=id.centers)

    kdat$assigned = kclust$clustering
    #---------------------------------------------------------------------------
    #=step 2: Compute the distance matrix between each point and the medoids
    #---------------------------------------------------------------------------

    # compute mean distance between centres
    distMean <- mean(stats::dist(centers))

    # transform back to lon/lat
    centers %<>%
      as.data.frame(row.names=NULL) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
      sf::st_transform(crs=sf::st_crs(4326)) %>%
      sf::st_coordinates()

    # calculate the distance to the center
    for(j in 1:nrow(centers)){
      kdat <- kdat %<>%
        dplyr::mutate(kdist = geosphere::distGeo(cbind(lon,lat), centers[j,]))
      # rename the column
      kdat <- kdat %>%
        dplyr::rename(!!paste0('kdist',j) := kdist)
    }
    if(iter > 0){
      # assign the outliers from each group to the closest cluster center
      switch(method[1],

             "top.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                   )
               }
             },

             "iqr.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                  .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                            .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             },

             "utmzone.outliers" = {
               j = 1
               kdat.copy <- kdat %>%
                 dplyr::filter(assigned == j) %>%
                 dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                 dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                 dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
               for(j in 2:k){
                 kdat.copy %<>%
                   rbind(kdat %>%
                           dplyr::filter(assigned == j) %>%
                           dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                           dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                           dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                   )
               }
             }
      )
      kdat.copy %<>% data.frame(row.names=NULL)
      kdat <- kdat.copy

    }
    #---------------------------------------------------------------------------
    #=step 3: Re-assign points to clusters
    #---------------------------------------------------------------------------
    kdat$index = 1:nrow(kdat)
    working = kdat
    ReAssignment = nrow(kdat) - (nrow(kdat) %% k)

    for(i in 1:ReAssignment){
      #cluster counts can be off by 1 due to uneven multiples of k.
      j = if(i %% k == 0) k else (i %% k)
      itemloc =  working$index[which.min(working[,(paste0("kdist", j))])][1]
      kdat$assigned[kdat$index == itemloc] = j
      working %<>%
        dplyr::filter(!index == itemloc)
    }
    # if there are some leftover points, assign to whoever's closest, without regard to k
    if(length(working$index)>0L){
      for(i in working$index){
        kdat$assigned[kdat$index == i] = which.min(working[working$index == i,grepl("kdist",colnames(working))])
      }
    }

    # assign the outliers from each group to the closest cluster center
    switch(method[1],

           "top.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::arrange_at(paste0("kdist", j), dplyr::desc) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(dplyr::row_number() <= n.outliers, is_kdist_min, assigned))
                 )
             }
           },

           "iqr.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate_at(.vars= dplyr::vars(!!paste0("kdist", j)),
                                          .funs= list(is.outliers = ~ . > stats::quantile(., probs=.75) + 1.5 * stats::IQR(.))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           },

           "utmzone.outliers" = {
             j = 1
             kdat.copy <- kdat %>%
               dplyr::filter(assigned == j) %>%
               dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
               dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
               dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
             for(j in 2:k){
               kdat.copy %<>%
                 rbind(kdat %>%
                         dplyr::filter(assigned == j) %>%
                         dplyr::mutate(is.outliers = LL2UTMZone(lon,lat,add_latitude_band=TRUE) != names(modal(LL2UTMZone(lon,lat,add_latitude_band=TRUE)))) %>%
                         dplyr::mutate(is_kdist_min =  dplyr::select(., !!paste0("kdist", 1:k)) %>% min.col()) %>%
                         dplyr::mutate(assigned = ifelse(is.outliers, is_kdist_min, assigned))
                 )
             }
           }
    )
    kdat.copy %<>% data.frame(row.names=NULL)
    kdat <- kdat.copy
    if(verbose){
      cat('#---------------------------------------------------------------------------\n')
      cat('#=step 4: Recalculate the medoids\n')
      cat('#---------------------------------------------------------------------------\n')
    }
    id.centers <- c()
    j = 1
    NewCenters <- kdat %>%
      dplyr::filter(assigned == j) %>%
      dplyr::select(lon, lat) %>%
      sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
      sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
      sf::st_coordinates() %>%
      cluster::pam(1) %$% list(medoids,id.med)

    id.centers <- c(id.centers, NewCenters[[2]])
    NewCenters <- NewCenters[[1]]

    while(j < k){
      j = j+1
      pamCenters <- kdat %>%
        dplyr::filter(assigned == j) %>%
        dplyr::select(lon, lat) %>%
        sf::st_as_sf(coords=1:2, crs=st_crs(4326)) %>%
        sf::st_transform(crs=sf::st_crs("+proj=robin")) %>%
        sf::st_coordinates() %>%
        cluster::pam(1) %$% list(medoids,id.med)

      id.centers <- c(id.centers, pamCenters[[2]])
      NewCenters %<>% rbind(pamCenters[[1]])
    }
    # New centroids
    NewCenters %<>% data.frame(row.names=NULL)
    centers <- NewCenters

    # keep assigments lon, lat and assigments only
    kdat <- kdat %>%
      dplyr::select(lon, lat, assigned)

    if(plot){
      cent <- centers %>%
        as.data.frame(row.names=NULL) %>%
        sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
        sf::st_transform(crs=sf::st_crs(4326)) %>%
        sf::st_coordinates() %>%
        as.data.frame(row.names=NULL)

      if(!is.factor(kdat$assigned))
        kdat$assigned %<>% as.factor()
      print(
        kdat %>% ggplot2::ggplot(ggplot2::aes(x = lon, y = lat, color = assigned)) +
          ggplot2::theme_minimal() + ggplot2::geom_point()  +
          ggplot2::geom_point(data = cent, ggplot2::aes(x = X, y = Y), color = "black", size = 4)
      )
    }

    process = abs(mean(stats::dist(centers)) - distMean) > tolerance.iter
    iter = iter + 1
  }
  # end iterations

  # assignment to factor
  if(!is.factor(kdat$assigned))
    kdat$assigned %<>% as.factor()

  centers %<>%
    sf::st_as_sf(coords=1:2, crs=st_crs("+proj=robin")) %>%
    sf::st_transform(crs=sf::st_crs(4326)) %>%
    sf::st_coordinates()

  return(list(Data=kdat, Centers=centers, converged=!process))
}

