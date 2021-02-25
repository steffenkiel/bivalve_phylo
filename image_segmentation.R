
library(stringr)
library(imager)
library(magick)


#--------------------------------------------------------------------------
# segments an image and saves the segments
# INPUT:
#  im_name <- path to a JPG-image
#  magick <- if TRUE, images are opened with magick's image_read function,
#            which handles large image files (> 3Mb) better than imager's 
#            load.image function. However, sometimes magick would not open
#            a file (no idea why), which never happened with imager.
# OUTPUT:
#  the segments saved in the same directory as the original image, 
#  with a running three-digit number added to the original file name
#--------------------------------------------------------------------------
segment_and_save <- function(im_name, magick=F)
{
  d <- dirname(im_name)
  
  # load image with magick
  if (magick){
    im <- image_read(im_name)
    if (image_info(im)$filesize > 1000000){
      im <- image_resize(im, geometry_size_pixels(width = 1000, height = 1000, preserve_aspect = TRUE))
    }
    im <- image_normalize(im)
    im <- magick2cimg(im)
  } else {
    # load image with imager
    im <- load.image(im_name)
  }
  
  # get coordinates of items
  im_h <- height(im)
  im_w <- width(im)
  xywh <- get.coordinates(im)
  
  if (length(xywh)>1){
    n_items <- nrow(xywh)          # number of items on the image
    name_base <- gsub(".jpg$|.JPG$", "_",im_name)
    new_names <- vector(length = n_items)
    
    for (i in 1:n_items){
      # get absolute coordinates for item on source image[i]
      absco <- absolute.coordinates(im_w,im_h,xywh[i,1],xywh[i,2],xywh[i,3],xywh[i,4])
      
      # copy specimen from original image
      item_im <- imsub(im, x>absco[1] & x<(absco[1]+absco[3]), y>absco[2] & y<(absco[2]+absco[4]))
      
      # create image name and save image [image name includes FULL PATH]
      new_names[i] <- paste(name_base, sprintf("%03d",i), ".jpg", sep = "")
      save.image(item_im,new_names[i],quality = 1)
    }
  } else {
    new_names <- NA
  }
  
  return(new_names)
}


#___________________________________________
#            helper functions

#-----------------------------------------------------------------
# segments the image and returns the coordinates of the segments
#-----------------------------------------------------------------
get.coordinates <- function(img)
{
  # thresholding
  wi <- grayscale(resize_fixed(img, 500))
  tim <- threshold(wi, adjust = .5)             # thresholding to identify items
  sp <- split_connected(tim) %>% purrr::discard(~ sum(.) < 200) # isolate & get rid of small stuff
  n.items <- length(sp)   # number of items
  
  if (n.items > 0){
    # get height and width of re-scaled image
    im_h <- height(wi)
    im_w <- width(wi)
    
    # make bounding boxes around each item and enlarge them a little
    items <- lapply(sp, bbox)
    items <- lapply(items, grow, 5)
    
    # make df for item data
    img.df <- as.data.frame(matrix(data = NA, nrow = n.items, ncol = 4))
    colnames(img.df) <- c("rel_x","rel_y","rel_w", "rel_h")
    
    # save coordinates for each bounding box
    for (i in 1:n.items){
      item <- where(items[[i]])
      img.df$rel_x[i] <- item[1,1] / im_w
      img.df$rel_y[i] <- item[1,2] / im_h
      img.df$rel_w[i] <- (item[nrow(item),1] - item[1,1]) / im_w
      img.df$rel_h[i] <- (item[nrow(item),2] - item[1,2]) / im_h
    }
  } else {
    img.df <- NA
  }
  
  return(img.df)
}


#--------------------------------------------------------
# turns coordinates relative to image size into
# x1, y1, w, h
#--------------------------------------------------------
absolute.coordinates <- function(im_w, im_h, item_x1, item_y1, item_w, item_h)
{
  out_x1 <- item_x1 * im_w
  out_y1 <- item_y1 * im_h
  out_width <- item_w * im_w
  out_height <- item_h * im_h
  
  out <- c(out_x1,out_y1,out_width,out_height)
  return(out)
}


#-------------------------------------------------------------------------
# converts an "imager" image to a "magick" image
#-------------------------------------------------------------------------
cimg2magick <- function(im, format = "png")
{
  tmp <- tempfile(fileext = paste(".", format, sep = ""))
  on.exit(unlink(tmp))
  imager::save.image(im, tmp)
  magick::image_read(tmp)
}
