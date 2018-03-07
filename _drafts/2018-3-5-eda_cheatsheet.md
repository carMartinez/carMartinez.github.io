---
layout: post
title: "EDA Cheatsheet - Python vs R"
permalink: /eda-cheatsheet/
---

A reference guide for exploratory data analysis (EDA) in Python and R. We'll be
using  genetic data from The Cancer Genome Atlas (TCGA) to walk through some routine
tasks when you first get your hands on a new (semi-clean) dataset.

# Motivation

R and Python both have their pros and cons. Hopefully you're able to pick the right
tool for the job. Sometimes I find myself using one language exclusively for
an extended period of time, and switching back to the other often entails some
re-learning of basic usage. This cheatsheet will give side-by-side comparisons of
R and Python code to jump-start the relearning process or act as a quick
reference guide for those unfamiliar with one or both languages.  

For the comparisons, I'll post code snippets below side-by-side
with Python on the left and R on the right. If you're interested in only a single
language, check out the Jupyter notebook and R markdown files for the original
Python and R code, respectively.

# The Dataset: TCGA

TCGA was a monumental, multi-institute research endeavor to characterize the genetic
landscapes of a variety of cancers. Essentially, hundreds of patients had tumor samples
profiled for genetic mutations. Most of the data is published online, and we'll be using
some publicly available data for this project. You don't have to be too familiar in
genetics to follow along with this analysis. We'll just explore some basic details
like which genes had the most mutations and how many mutations existed per tumor.
I'll give some simple explanations where necessary.    

The data itself can be downloaded from:

I will not the host the data myself, but I can provide the manifest file needed
to download the same data used here. Specifically, I used all publicly
available data from the breast cancer (BRCA) dataset. However, the analysis
is quite general so you could use any dataset from TCGA that is in the `.maf`
format (detailed here & here). We'll see below where to define the input paths
to the data, so you can change the paths to suit whichever dataset you download (or
adapt it to an unrelated dataset).


# Dependencies

First, let's load the dependencies. For Python we'll be making heavy use of the
`pandas` library for data wrangling plus `matplotlib` and `seaborn` (built on
top of `matplotlib`) for visualizations. Most of the R packages needed are
contained in the `tidyverse` mega-package.

<div class="two_col_wrap">
<div class="left_code">
{% highlight python %}
import numpy as np
import pandas as pd
from collections import OrderedDict
from IPython.display import display

import matplotlib.pyplot as plt
import seaborn as sns
%matplotlib inline
{% endhighlight %}

</div>
<div class="right_code">

{% highlight R %}
library(tidyverse)
library(magrittr)
{% endhighlight %}

</div>
</div>
<div style="clear: both;"></div>


# Reading Data

Both `pandas` and `readr` have the capability to read `gzip` files which is
convenient as the downloaded `.maf` files are compressed by default. It's worth
noting that the `readr` package is significantly faster at reading in files than
base R, so use it whenever possible (or `fread` from `data.table`, but we won't
cover that here).

# Text Summaries
