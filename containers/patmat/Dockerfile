FROM ubuntu:22.04
ARG DEBIAN_FRONTEND=noninteractive

LABEL \
  version="${PATMAT_VERSION}" \
  description="Patmat image to use in the POAGA pipeline"

ENV SAMTOOLS_VERSION=1.21
ENV PATMAT_VERSION=dev

RUN apt-get update
RUN apt-get install -y --no-install-recommends software-properties-common dirmngr gdebi-core python3-pip python3-dev git bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev libssl-dev libcurl4-openssl-dev autoconf automake \
  && cd /usr/local/bin \
  && ln -s /usr/bin/python3 python \
  && pip3 --no-cache-dir install --upgrade pip

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd samtools-${SAMTOOLS_VERSION}/ \
	&& ./configure && make && make install

WORKDIR /opt
RUN wget https://github.com/samtools/bcftools/releases/download/${SAMTOOLS_VERSION}/bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf bcftools-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf bcftools-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd bcftools-${SAMTOOLS_VERSION}/ \
	&& ./configure && make && make install

WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/${SAMTOOLS_VERSION}/htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& tar -xjf htslib-${SAMTOOLS_VERSION}.tar.bz2 \
	&& rm -rf htslib-${SAMTOOLS_VERSION}.tar.bz2  \
	&& cd htslib-${SAMTOOLS_VERSION}/ \
	&& ./configure && make && make install

WORKDIR /opt
RUN git clone --recursive https://github.com/epi2me-labs/modbam2bed.git && cd modbam2bed && make

# Bedtools
RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz \
        && tar -zxvf bedtools-2.29.1.tar.gz \
        && cd bedtools2 \
        && make && make install

# Modkit
RUN wget https://github.com/nanoporetech/modkit/releases/download/v0.4.1/modkit_v0.4.1_u16_x86_64.tar.gz \
        && tar xf modkit_v0.4.1_u16_x86_64.tar.gz \
        && cp dist_modkit_v0.4.1_cec0a0b/modkit /usr/local/bin

# R
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc
RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" -y
RUN apt install -y r-base r-base-core r-recommended r-base-dev
RUN add-apt-repository ppa:c2d4u.team/c2d4u4.0+ -y
RUN apt install -y r-cran-ggplot2 r-cran-biocmanager r-cran-tibble
RUN apt install -y r-bioc-beachmat
RUN Rscript -e 'install.packages(c("argparse", "sys", "R.utils"), repos="http://cran.us.r-project.org")'

RUN pip3 install pysam pytabix tqdm modbampy

# Patmat
WORKDIR /opt
RUN git clone https://github.com/vahidAK/PatMat && cd PatMat \
	&& git checkout $PATMAT_VERSION \
    && ln -s $PWD/{patmat,Strand-seq}/* /opt

# Now we need to open R and install the required R packages as follows:
RUN Rscript -e 'install.packages(c("argparse", "sys", "R.utils"), repos="http://cran.us.r-project.org")'
RUN Rscript -e 'install.packages(c("pkgconfig", "pkgdepends", "pak", "curl"), repos="http://cran.us.r-project.org")'
RUN Rscript -e 'library(pkgdepends); library(pkgconfig); library(pak); pak::pkg_install("vincent-hanlon/InvertypeR")'
RUN Rscript -e 'BiocManager::install("DSS")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'

# Some cleanup
RUN rm -rf /var/lib/apt/lists/*




