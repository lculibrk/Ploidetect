## Based on r-base dockerfile
FROM rocker/r-base:3.6.3

## Install dependencies
RUN apt-mark hold r-base-core
RUN  apt-get update -qq && apt-get -y --no-install-recommends install \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libglib2.0-dev \
  libmount-dev \
  libblkid-dev \
  libgtk2.0-dev \
  libgdk-pixbuf-2.0-dev \
  libpango1.0-dev \
  libtiff-dev \
  libharfbuzz-dev \
  libcairo2-dev \
  libcairo2-dev \
  xvfb \
  xauth \
  xfonts-base \
  libfontconfig1-dev \
  libfreetype6-dev \
  libfreetype-dev \
  libfontconfig-dev \
  build-essential gcc g++ cpp g++-10 gcc-10-base cpp-10 \
  libxt-dev

## Install devtools & Ploidetect
RUN install2.r devtools tidyr ggrastr 
RUN Rscript -e "devtools::install_github('lculibrk/Ploidetect', ref = 'v1.4.2')"

## Install preprocessing software
#RUN apt-get -y --no-install-recommends install \ 
#  samtools \ 
#  bedtools \ 
#  bcftools \ 
#  python3-pip
SHELL ["/bin/bash", "-c"]
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh
RUN bash Mambaforge-$(uname)-$(uname -m).sh -b -p /mamba/
#RUN source /mamba/bin/activate
ENV PATH=/mamba/bin:${PATH}
RUN mamba init bash
#RUN echo "source /mamba/bin/activate" >> ~/.bashrc
RUN mamba install -yc bioconda samtools=1.13 bcftools bedtools


RUN Rscript -e "install.packages('BiocManager')"
RUN Rscript -e "BiocManager::install('GenomicRanges')"
RUN chmod +x /usr/local/lib/R/site-library/Ploidetect/plot_regions.R
ENV PATH="/usr/local/lib/R/site-library/Ploidetect/:${PATH}"

ENTRYPOINT ["bash"]
## pip install the required dependencies
#RUN pip3 install fileinput