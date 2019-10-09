# get the base image, the rocker/verse has R, RStudio and pandoc
FROM rocker/verse:3.5.2

# required
MAINTAINER Hao Ye <lhopitalified@gmail.com>

COPY . /portalDS

# go into the repo directory
RUN apt-get update \
  && sudo apt-get install -y --no-install-recommends \
    libudunits2-dev \
    libgsl0-dev \
  # build this compendium package
  && R -e "devtools::install('/portalDS', dep=TRUE)"
