---
layout: post
title: "NAIP Quadrangle Downloader"
permalink: /naip-quadrangle-downloader/
description: "Handy script to automate downloading of NAIP aerial images from AWS based only
on simple names like 'Statue of Liberty, NY'"
---

![Four US Landmarks](/assets/img/naip_downloader/four_us_landmarks.png)

The National Agriculture Imagery Program (NAIP) provides high-resolution aerial imagery of the United States to the public via an [Amazon Web Services (AWS) S3 bucket](https://aws.amazon.com/public-datasets/naip/).

Essentially, the US is broken into a grid of cells called quadrangles,
and the images NAIP provides correspond to cells in this grid. Thus,
to get an image containing a specific city or landmark, one must first
figure out in which quadrangle it is located. This can be a little tedious
to do by hand, so I made a function `download_landmark()` to automate it.

{% highlight python %}
download_landmark('Statue of Liberty, NY')
{% endhighlight %}

This handy function can easily be incorporated into geographic information systems
workflows. For example, although the image downloaded indeed contains the Statue
of Liberty, it also has a lot more geography and you'd have to go searching through it
to find the landmark. With only a bit more work, we can use the function not
just to download the image but also extract a subset of it centered around the
landmark. In fact, we'll do it for a few landmarks in the code below which
produces the image at the top of this page.

{% highlight python %}
import rasterio
from rasterio.mask import mask as rio_mask
from matplotlib import pyplot as plt
from shapely.geometry import box
from pyproj import Proj
import numpy as np

places = [
    'Statue of Liberty, NY',
    'National Mall, Washington DC',
    'Hoover Dam, AZ',
    'Golden Gate Bridge, CA'
]
img_dir = 'naip/img'
mask_len = 1500 # meters
fig, axes = plt.subplots(2, 2, figsize = (10, 10))
axes = axes.flatten()

for i, place in enumerate(places):
    img_path, lat, lon = download_landmark(place, img_dir)
    with rasterio.open(img_path) as src:

        ## Project landmark lat,lon
        proj = Proj(src.crs)
        x,y = proj(lon, lat)

        ## Make mask centered at the landmark
        mask_shape = box(x-mask_len/2, y-mask_len/2, x+mask_len/2, y+mask_len/2)
        mask = {'type': 'Polygon', 'coordinates': [list(mask_shape.exterior.coords)]}

        ## Extract only the data within the mask
        img, transform = rio_mask(src, [mask], crop = True)
        img = np.moveaxis(img, 0, 2) # channels last

        axes[i].set_title(place)
        axes[i].imshow(img)
fig.savefig('four_us_landmarks.png', bbox_inches = 'tight')
{% endhighlight %}

The script as well as a tutorial explaining it can be found on my
[GitHub](https://github.com/carMartinez/naip_quadrangle_downloader).
