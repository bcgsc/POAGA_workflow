FROM condaforge/mambaforge:latest

LABEL \
  version="${PATMAT_VERSION}" \
  description="Patmat image to use in the POAGA pipeline"

ENV SAMTOOLS_VERSION=1.21
ENV PATMAT_VERSION=v1.4.0

ENV TZ=America/Vancouver
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone
RUN apt-get update && apt-get install -y  bzip2 gcc g++ make zlib1g-dev wget libncurses5-dev liblzma-dev libbz2-dev libssl-dev libcurl4-openssl-dev autoconf automake && rm -rf /var/lib/apt/lists/*

WORKDIR /opt
RUN wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        && tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2 \
        && rm -rf samtools-${SAMTOOLS_VERSION}.tar.bz2  \
        && cd samtools-${SAMTOOLS_VERSION}/ \
        && ./configure --without-curses && make && make install

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
RUN apt update -qq \
        && apt install -y  --no-install-recommends software-properties-common dirmngr \
        && wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | tee -a /etc/apt/trusted.gpg.d/cran_ubuntu_key.asc \
        && add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" \
        && apt install -y --no-install-recommends r-base


WORKDIR /opt
RUN git clone https://github.com/vahidAK/PatMat && cd PatMat \
	&& git checkout $PATMAT_VERSION \
    && ln -s $PWD/{patmat,Strand-seq}/* /opt

RUN conda init && cd /opt/PatMat && chmod +x Strand-seq/*.R && conda env create -f patmat_env.yml
SHELL ["conda", "run", "-n", "patmat", "/bin/bash", "-c"]
RUN  cd /opt/PatMat && poetry install

RUN which R && R --version
RUN Rscript -e 'install.packages(c("argparse", "sys", "R.utils"), repos="https://cloud.r-project.org")'
RUN Rscript -e 'install.packages(c("ggplot2", "BiocManager", "tibble"), repos="https://cloud.r-project.org")' && Rscript -e 'BiocManager::install("beachmat")'

#RUN Rscript -e 'install.packages(c("pkgconfig", "pkgdepends", "pak", "curl"), repos="http://cloud.r-project.org")' && Rscript -e 'library(pkgdepends); library(pkgconfig); library(pak); pak::pkg_install("vincent-hanlon/InvertypeR")'
RUN Rscript -e 'remotes::install_github("vincent-hanlon/InvertypeR")'

RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")'
RUN Rscript -e 'BiocManager::install("BSgenome.Hsapiens.UCSC.hs1")'
RUN Rscript -e 'BiocManager::install("DSS")'

ENV PATH=/usr/local/bin/:/opt/conda/envs/patmat/bin:$PATH
