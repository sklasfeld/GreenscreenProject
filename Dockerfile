FROM ubuntu:latest

MAINTAINER Samantha Klasfeld <sjk314@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update --fix-missing && \
	apt-get install -y --no-install-recommends \
	build-essential zip wget git \
	default-jre bedtools \
	r-base r-cran-randomforest \
	python3.6 python3-pip python3-setuptools python3-dev \
	vim nano less \
	libncurses5-dev zlib1g-dev libbz2-dev \
	liblzma-dev libcurl3-dev libxml2-dev && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
	
# destination to install specific softwares
ENV DEST=/usr/src

WORKDIR $DEST

# install SRA toolkit
RUN wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
	tar -xzvf sratoolkit.current-ubuntu64.tar.gz

ENV PATH=${PATH}:/usr/src/sratoolkit.2.11.1-ubuntu64/bin
# install FASTQC

RUN wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip -O FastQC-0.11.9.zip && \
	unzip FastQC-0.11.9.zip

RUN chmod +x /usr/src/FastQC/fastqc

ENV PATH=${PATH}:/usr/src/FastQC

# install trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip -O Trimmomatic-0.39.zip && \
	unzip Trimmomatic-0.39.zip 

# install bowtie2
RUN wget https://github.com/BenLangmead/bowtie2/releases/download/v2.4.4/bowtie2-2.4.4-linux-x86_64.zip -O bowtie2-2.4.4-linux-x86_64.zip && \
	unzip bowtie2-2.4.4-linux-x86_64.zip

ENV PATH=${PATH}:/usr/src/bowtie2-2.4.4-linux-x86_64

# install samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
	tar jxf samtools-1.9.tar.bz2 && \
	rm samtools-1.9.tar.bz2 && \
	cd samtools-1.9 && \
	./configure --prefix $(pwd) && \
	make

ENV PATH=${PATH}:/usr/src/samtools-1.9

# install PICARD
RUN git clone https://github.com/broadinstitute/picard.git

WORKDIR $DEST/picard

RUN $DEST/picard/gradlew shadowJar

# reset working directory

WORKDIR /

# install python libraries

COPY requirements.txt ./

RUN pip install --no-cache-dir --upgrade pip && \
    pip install --no-cache-dir -r requirements.txt

# install R libraries

RUN Rscript -e "install.packages('argparse')"

RUN Rscript -e "install.packages('BiocManager')"

RUN Rscript -e "BiocManager::install('ChIPQC')"

RUN Rscript -e "BiocManager::install('GenomicFeatures')"

RUN Rscript -e "BiocManager::install('GenomicRanges')"

# copy scripts and helper meta information into container

COPY scripts /home/scripts

COPY meta /home/meta

# report done

CMD ["echo", "docker", "image", "built"]

WORKDIR /home