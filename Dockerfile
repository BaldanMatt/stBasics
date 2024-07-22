# base image
FROM rstudio/r-base:4.3-jammy

# set working directory in container
# WORKDIR /app

# copy the current directory contents into the container at /app
# COPY . /app

# Install needed system libraries
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
		libudunits2-dev \
		libgdal-dev \
		libgeos-dev \
		libproj-dev \
		libfftw3-dev \
		libmagick++-dev

# Install R packages
RUN R -e "install.packages('BiocManager')"
RUN R -e "BiocManager::install('SpatialFeatureExperiment', update=TRUE, ask=FALSE)"

# Make port available to run R
EXPOSE 8787
