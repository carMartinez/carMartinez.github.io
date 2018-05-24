---
layout: post
title: "Perimeter Extraction"
permalink: /perimeter-extraction/
---

Automatically extract the perimeter of a city given an aerial or satellite image
of the city. A convolutional neural network is used to identify the road network
in the image, and then image processing techniques expand the road network into
a contour of the city.

![Method diagram](/assets/img/perimeter/summary.png)

{:.image-caption}
Figure 1. Diagram of extraction method. a) Aerial image containing Lake Havasu
City, AZ. b) Road network prediction from a CNN. c) Mesh surrounding the city after
morphological operations to expand the road network. d) Contour of the largest
mesh object overlaid on the original aerial image.

# Overview

The process for perimeter extraction (see Figure 1) is largely inspired from
[this paper](https://ieeexplore.ieee.org/document/7566791/).
In it, the authors construct Google Map-esque images of cities using freely
available Open Street Map (OSM) data. The road networks in these images are detected with
elementary image processing techniques, and morphological operations expand
the road networks into binary masks that encompasses the cities. A simple contour procedure
around these masks produce the perimeters.

In this project, we accomplish the same perimeter detection task but
with aerial images instead of artificial Google Maps-esque images. The overall scheme
is essentially the same: identify the road network and apply morphological and
contour operations to create a bounding perimeter. Indeed, the morphological and
contour operations shown in Figure 1 are essentially identical in the reference paper.
The significant difference comes in the road detection step.

In the reference paper, road detection was relatively simple as the authors were
able to choose specific colors for roads when constructing the synthetic images. Thus,
road detection could be easily accomplished by filtering for colors or even with
edge detectors
operating on relevant color channels. Such techniques are not sufficient for road
detection in aerial images as the source is not so nice: roads can vary in color
and may blend in with their surroundings, and edge detectors will pick up noise
such as buildings or natural scenery. Thus, a more sophisticated method is necessary.

Other work has shown that machine learning can be successfully used for road
detection and the more general task of image segmentation. Here we use a convolutional
neural network (CNN) with a modified U-Net architecture to perform road detection
on aerial images. Training and
benchmark data is generated via OSM road labels.

The results of the procedure for 5 Arizona cities (Bullhead, Flagstaff, Globe,
Lake Havasu City, and Yuma) are presented in the [Results section](#results). The
outline for data acquisition, CNN training, and perimeter extraction procedure are
detailed in the [Methods section](#methods).

# Results

# Methods

This section outlines the procedure for this project from data acquisition to
CNN training to application. A more thorough, low-level
tutorial of each step including full code is given in individual Jupyter
notebooks on GitHub (one notebook per subsection below).

All computations were performed on a desktop computer with an Intel i5 6600K processor,
16GB of RAM, and an Nvidia 1070ti GPU.

## Dataset

The georeferenced aerial images are sourced from the US National Agriculture Imagery Program (NAIP) which
generated 1-meter-per-pixel aerial images of the entire US region. The data is freely
available on an Amazon Web Services (AWS) S3 bucket (bandwidth fees apply if you are not in the US-East region).

For this project we focus on cities in Arizona. In particular, aerial images for
Phoenix, Bullhead, Lake Havasu, Flagstaff, Globe, and Yuma are downloaded.

A single NAIP image covers only a small portion of land (based on the USGS quadrangle specification), so for each
city multiple images must be downloaded. The Google Maps' geolocator service is used to determine
an appropriate bounding box for each city, and then NAIP images
are downloaded to cover the entire bounding box.

## Labels

Freely available OSM data is used to label the roads in all the images downloaded
in the previous section. OSM provides shapefiles which contain line geometries for
roads throughout the world. The latitude, longitude coordinates of the roads are
matched to the georeferenced aerial images to produce binary mask images where pixels
are labeled 1 if they contain a road and 0 otherwise.

However, OSM road geometries correspond to road centers rather than the full width
of the roads. To compensate, road labels are expanded with a 3x3 kernel dilation
operation applied twice (see Figure 2).

<table class="table-center"><tr><td><img class="img-responsive" src='/assets/img/perimeter/labels_no_dilation.png'></td><td><img class="img-responsive" src='/assets/img/perimeter/labels_dilation.png'></td></tr></table>

{:.image-caption}
Figure 2. Left: Road labels without dilation. Right: Twice dilated.

Although road labels are produced for all images, only data corresponding to the
Phoenix, AZ area are used for training and evaluating the CNN performance.
Since the Phoenix area is rather large and neither images nor labels can
comfortably fit in 16GB of memory, the Phoenix data is split into Training/Development/Testing
partitions and exported to an HDF5 dataset.

## Preprocess

For implementing the perimeter detection procedure, it is convenient to have a single
image that encompasses a city rather than multiple smaller images. In this step,
the GDAL program is used to mosaic the multiple images previously downloaded into
a single georeferenced image for each city (as in Figure 1a).

## Learn

The aerial images and respective labels corresponding to the Phoenix, AZ area were
previously organized into training/development/test sets in the [Labels](#labels) section.

A modified U-Net CNN architecture was fit to the training set with performance metrics
evaluated on the development set for refinement. The CNN takes as input a 640x640x3-pixel RGB image and produces a 512x512x1 image with pixel values corresponding to
the probability of that pixel containing a road. The architecture here is very similar
to the original U-Net but with added dropout layers and summation skip connections rather than concatenation. The
layer-by-layer summary is reported in the Jupyter notebook. Final metrics were
evaluated on the test set.

It is worth noting that the U-Net outputs a different sized image than what is input. To
account for this, a large RGB image is split into 512x512-pixel tiles. All these tiles are padded with 64-pixel reflections
on each side to produce the 640x640 input. The final 512x512 output corresponds to the
original 512x512 input dimensions.

The model was implemented via the Keras library with Tensorflow GPU backend.

## Predict

The CNN trained in the previous section was applied to five Arizona cities:
Bullhead, Lake Havasu, Flagstaff, Globe, and Yuma. For each city, the mosaic
created in the preprocessing step was split into 512x512x3 tiles each of which
were then padded
with 64-pixel reflections on every side to form the 640x640x3 inputs needed. The
resulting prediction tiles were rearranged back into a single image roughly the
same size as the original mosaic (see Figure 1b); some pixels were cropped off the mosaic to
ensure it could be evenly divided into the 512x512-pixel tiles.   

## Image Processing

The prediction images in the last step were rounded to produce a single binary
mask for each city with pixel values of 1 indicating the presence of a road
and 0 otherwise. Morphological operations via the OpenCV library were applied to
these images.

First, a small dilation with a square 3x3 kernel was applied to remove small
prediction artifacts. Then, the road network was expanded through a dilation
operation repeated 15 times with a large 15x15 kernel. The goal of this operation
is to connect unconnected roads and make the road network into a single,
completely connected mesh that surrounds the entire city. Lastly, a dilation
with the same 15x15 kernel was applied and repeated 15 times. This ensures the
boundary of the road network mesh does not exceed the bounds of the original,
undilated predictions. See Figure 1c for the resulting mesh.

After the morphological operations were applied, contours were found for all
connected objects in the image (as in Figure 1c). Often there were multiple
distinct, unconnected
objects, so the contour of the largest object was taken to be the city perimeter.
This contour could then be overlaid on the original aerial image for visual
assessment (Figure 1d).
