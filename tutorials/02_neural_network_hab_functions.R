#' Takes raw input data, filters for image criteria and creates images with dimensions (n_steps + forecast steps) x length(toxins + environmentals)
#' 
#' @param raw_data database with toxin measurements with their date sampled, location, shellfish species and additional environmental data
#' @param tox_levels toxin level categories used for classifying total toxicity
#' @param forecast_steps the number of weeks ahead of the image the forecast is made for
#' @param n_steps the number of weeks of samples in an image
#' @param minimum_gap the smallest gap between samples allowed into an image
#' @param maximum_gap the largest gap between samples allowed into an image
#' @param toxins list of individual paralytic toxin names (12) for toxin columns
#' @param environmentals environmental variables
#' @return each list is an image along with its associated data (location_id, date, etc.)
#' \itemize{
#' \item{status logical if the image passes the image gap criteria (gap >= minimum_gap & gap <= maximum_gap)}
#' \item{year the year the image is from}
#' \item{location_id the sampling station id}
#' \item{toxixty the total toxicity used to regress on instead of classify a binned toxicity}
#' \item{classification the classification (0:num_classes) of the final row in the image}
#' \item{date the date of the final row in the image (the forecast is for forecast_steps ahead of this date)}
#' \item{image a 2 dimensional array with the dimensions (n_steps + forecast steps) x length(toxins + environmentals)}
#' }
#' 
#' @export
make_image_list <- function(raw_data, tox_levels, forecast_steps, n_steps, minimum_gap, maximum_gap, toxins, environmentals) {
  
  stations <- pspdata::read_station_metadata()
  
  lon_lut <- stations$lon
  names(lon_lut) <- stations$location_id
  
  normalized_data <- raw_data %>% 
    dplyr::mutate(classification = recode_classification(.data$total_toxicity, tox_levels),
                  meets_gap = check_gap(.data$gap_days, minimum_gap, maximum_gap),
                  week = date_to_week(.data$date),
                  lon = lon_lut[.data$location_id]) %>% 
    normalize_input(toxins, environmentals)
  
  find_images <- function(tbl, key, forecast_steps, n_steps, minimum_gap, maximum_gap, toxins, environmentals) {
    
    make_images <- function(batch, tbl, forecast_steps, n_steps, minimum_gap, maximum_gap, toxins, environmentals) {
      
      image_batch <- tbl %>% dplyr::slice(batch)
      
      if (any(image_batch$meets_gap[2:(n_steps+forecast_steps)] == FALSE)) {
        z <- list(status=FALSE)
      } else {
        image <- as.matrix(dplyr::ungroup(image_batch) %>% 
                             dplyr::select(dplyr::all_of(c(toxins, environmentals))))
        
        if (image_batch$id[n_steps+forecast_steps] == "FORECAST_WEEK") {
          z <- list(status=          TRUE,
                    year =           "FORECAST_IMAGE",
                    location_id =    image_batch$location_id[1],
                    classification = image_batch$classification[n_steps+forecast_steps],
                    toxicity =       image_batch$total_toxicity[n_steps+forecast_steps],
                    date =           image_batch$date[n_steps],
                    image =          image[1:n_steps,])
        } else {
          z <- list(status=          TRUE,
                    year =           image_batch$year[1],
                    location_id =    image_batch$location_id[1],
                    classification = image_batch$classification[n_steps+forecast_steps],
                    toxicity =       image_batch$total_toxicity[n_steps+forecast_steps],
                    date =           image_batch$date[n_steps],
                    image =          image[1:n_steps,])
          
        }
        
        return(z)
      }
    }
    
    if (nrow(tbl) < (n_steps+forecast_steps)) {
      return(NULL)
    }
    
    nbatches <- n_batches(nrow(tbl), (n_steps+forecast_steps))
    batches <- compute_batches(nbatches, (n_steps+forecast_steps))
    
    xx <- lapply(batches, make_images, tbl, forecast_steps, n_steps, minimum_gap, maximum_gap, toxins, environmentals)
    gap_verified <- sapply(xx, function(x){return(x$status)})
    
    xx <- xx[gap_verified]
    
    return(xx)
  }
  
  image_list <- normalized_data %>%
    dplyr::group_by(.data$location_id, .data$year) %>%
    dplyr::arrange(date) %>% 
    dplyr::group_map(find_images, forecast_steps, n_steps, minimum_gap, maximum_gap, toxins, environmentals, .keep=TRUE) %>% 
    unlist(recursive = FALSE)
  
  return(image_list)
  
}



#' Crates image and labels for input into neural net
#' Image takes all images in psp_lst and stretches them into an array
#' Labels takes classifications and categorizes them for nn input
#' 
#' @param image_list_subset subset of image list from make_image_list() for either training or testing data
#' @param num_classes the number of toxicity classification categories 
#' @param missing_value value to replace na with
#' @param scaling_factors null if training data; training data scaling factors are passed to scale testing data
#' @param scaling selected method to scale input data
#' @param downsample logical indicating whether or not to balance the frequency of each class in the training images
#' @param upsample logical to call a function that balances class distribution by upsampling rare classes
#' @return list containing the formatted images and their labels as keras model input, and additional data
#' \itemize{
#' \item{labels}
#' \item{image a 2 dimensional array where each row is an image and the columns are toxins and environmentals from each week}
#' \item{classifications the classification of each image}
#' \item{locations the sampling station of each image}
#' \item{dates the date of the final week of each image}
#' \item{scaling_factors}
#' }
#' 
#' @export
pool_images_and_labels <- function(image_list_subset, 
                                   num_classes = 4, 
                                   missing_value = 0.5, 
                                   scaling_factors = NULL, 
                                   scaling = c("normalize", "input_scale")[2],
                                   downsample=FALSE,
                                   upsample=FALSE) {
  
  xx <- unlist(image_list_subset, recursive = FALSE) 
  
  if (upsample == TRUE) {
    xx <- xx %>% 
      upsample()
  }
  
  if (downsample == TRUE) {
    xx <- xx %>% 
      balance_classes() %>% 
      sample()
  } else{
    xx <- xx %>% 
      sample()
  }
  
  dim_image <- dim(xx[[1]]$image)
  
  images <- lapply(xx, function(x){return(x$image)})
  
  # Replace any NA values with specified missing value
  # @param x 
  # @param missing_value to replace na toxin levels with
  # @return x
  replace_na <- function(x, missing_value = -1) {
    x[is.na(x)] <- missing_value
    return(x)
  }
  
  image <- abind::abind(images, along = 3) %>% 
    aperm(c(3, 1, 2)) %>% 
    reticulate::array_reshape(c(length(xx), prod(dim_image)))
  
  image <- image %>% 
    replace_na(missing_value = missing_value)
  
  labels <- sapply(xx, function(x){return(x$classification)}) %>% 
    keras::to_categorical(num_classes = num_classes)
  
  classifications <- sapply(xx, function(x){return(x$classification)})
  attr(classifications, "names") <- NULL
  
  locations <- sapply(xx, function(x){return(x$location_id)})
  attr(locations, "names") <- NULL
  
  dates <- sapply(xx, function(x){return(x$date)})
  attr(dates, "names") <- NULL
  
  toxicity = sapply(xx, function(x){x$toxicity})
  attr(toxicity, "names") <- NULL
  
  r <- list(labels = labels, 
            image = image, 
            classifications = classifications,
            toxicity = toxicity,
            locations = locations,
            dates = dates,
            scaling_factors = scaling_factors)
  
  return(r)
  
} 
