FROM rocker/tidyverse

LABEL maintainer="ckaipf@posteo.de"

# Install bedtools from Debian package
RUN apt-get update && apt-get install -y \
      bedtools 