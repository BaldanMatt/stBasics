---
title: "bioConductor"
author: "matteo baldan"
format: html
editor: visual
---

## Importing packages

```{r}
library(SpatialExperiment)
```

```{r}
example("read10xVisium", echo=FALSE)
```

```{r}
spe
```

# Spatial coords

These are stored as a numeric matrix of spatial coordinates (x, y and enventually z). They are stored inside the \`int_colData\` and are accessible by the corresponding accessor

```{r}
head(spatialCoords(spe))
```

The names are accessible as well by the corresponding accessor

```{r}
spatialCoordsNames(spe)
```

# Img Data

These are stored inside the int_metadata field called imgData as a DataFrame.\
The DataFrame has the following model:

-   each row is an image for a given sample with unique identifier

-   the column specify

    -   sample_id:

    -   image_id:

    -   data: SpatialImage object

    -   scaleFactor: adjusts pixel positions of the original, full resolution image to pixel posistions in the image

```{r}
imgData(spe)
```

# SpatialImage Class

Images are stored within data of imgData as a list of SpatialImages. Where eahc one may be on of the following:

-   LoadedSpatialImage: fully realized into memory as a raster. @image

-   StoredSpatialImage: image stored in a local file and loaded into memory only on request. @path

-   RemoteSpatialImage: remotely hosted and retrieved through an URL on request. @url

```{r}
(spi <- getImg(spe))
identical(spi, imgData(spe)$data[[1]])
plot(imgRaster(spe))
```

```{r}
url <- "https://i.redd.it/3pw5uah7xo041.jpg"
spe <- addImg(spe, 
    sample_id = "section1", 
    image_id = "pomeranian",
    imageSource = url, 
    scaleFactor = NA_real_, 
    load = TRUE)
img <- imgRaster(spe, 
    sample_id = "section1", 
    image_id = "pomeranian")
plot(img)
```
