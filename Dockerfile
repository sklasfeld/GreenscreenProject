FROM ubuntu:latest

MAINTAINER Samantha Klasfeld <sjk314@gmail.com>

ENV DEBIAN_FRONTEND=noninteractive

#install dependencies 

RUN apt-get update --fix-missing && \
    apt-get install -y --no-install-recommends \
    bash-completion ca-certificates file \
    fonts-texgyre g++ gfortran gsfonts \
    libblas-dev libcurl3-dev libcurl4 \
    liblapack-dev liblzma5 liblzma-dev libncurses5-dev \
    libopenblas-dev libpangocairo-1.0-0 libpcre3 libpng16-16 \
    libssl-dev libssl-doc libtiff5 librsvg2-dev \
    libxml2-dev locales make tar unzip wget curl zip \
    zlib1g zlib1g-dev && \
    echo "en_US.UTF-8 UTF-8" >> /etc/locale.gen \
    && locale-gen en_US.utf8 \
    && /usr/sbin/update-locale LANG=en_US.UTF-8 \
    > /dev/null

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    default-jdk libbz2-dev \
    libcairo2-dev libcurl4-openssl-dev \
    libpango1.0-dev libjpeg-dev \
    libicu-dev libpcre3-dev \
    libpng-dev libreadline-dev libtiff5-dev \
    liblzma-dev libx11-dev librsvg2-bin libv8-dev \
    libxt-dev perl tcl8.6-dev tk8.6-dev \
    texinfo texlive-extra-utils texlive-fonts-recommended \
    texlive-fonts-extra texlive-latex-recommended \
    texlive-latex-base texlive-latex-extra \
    x11proto-core-dev xauth xfonts-base \
    xvfb zlib1g-dev >> /dev/null

# install R

WORKDIR /tmp

RUN wget https://cran.r-project.org/src/base/R-4/R-4.0.0.tar.gz

RUN tar -xf R-4.0.0.tar.gz

WORKDIR /tmp/R-4.0.0

## set R compiler flags
RUN R_PAPERSIZE=letter \
    R_BATCHSAVE="--no-save --no-restore" \
    R_BROWSER=xdg-open \
    PAGER=/usr/bin/pager \
    PERL=/usr/bin/perl \
    R_UNZIPCMD=/usr/bin/unzip \
    R_ZIPCMD=/usr/bin/zip \
    R_PRINTCMD=/usr/bin/lpr \
    LIBnn=lib \
    AWK=/usr/bin/awk \
    CFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g" \
    CXXFLAGS="-g -O2 -fstack-protector-strong -Wformat -Werror=format-security -Wdate-time -D_FORTIFY_SOURCE=2 -g"


RUN /tmp/R-4.0.0/configure --enable-R-shlib \
               --enable-memory-profiling \
               --with-readline \
               --with-blas \
               --with-tcltk \
               --disable-nls \
               --with-recommended-packages \
               --with-x=no \
    && make \
    && make install

RUN apt-get clean && apt-get purge && \
    rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
    >> /dev/null

# install R libraries
RUN R -e "install.packages('rsvg', dependencies=TRUE, repos='https://cran.r-project.org')"
RUN R -e "install.packages('BiocManager', dependencies=TRUE, repos='http://cran.us.r-project.org')"

RUN wget https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz && \
   tar -xzvf locfit_1.5-9.4.tar.gz
RUN R -e "install.packages('locfit', repos = NULL, type='source')"


RUN Rscript -e "BiocManager::install()"
RUN Rscript -e "BiocManager::install('ShortRead', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('edgeR', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('DESeq2', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('GOstats', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('amap', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('systemPipeR', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('ChIPQC', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('GenomicFeatures', force=TRUE, ask=FALSE)"
RUN Rscript -e "BiocManager::install('GenomicRanges', force=TRUE, ask=FALSE)"

# import python & other useful software
RUN apt-get update --fix-missing && \
	apt-get install -y --no-install-recommends \
	build-essential gzip git \
	default-jre bedtools \
	python3.6 python3-pip python3-setuptools python3-dev \
	vim nano less rsync && \
	apt-get clean && apt-get purge && \
	rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

# install python libraries

COPY requirements.txt ./

RUN pip3 install --no-cache-dir --upgrade pip && \
    pip3 install --no-cache-dir -r requirements.txt
	
# destination to install specific softwares
ENV DEST=/usr/src

WORKDIR $DEST

# install KentUtils

RUN mkdir kentUtils

WORKDIR $DEST/kentUtils

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedCoverage ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedExtendRanges ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedGraphToBigWig ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedIntersect ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedSort ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bedToBigBed ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigAverageOverBed ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigMerge ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigSummary ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigToBedGraph ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/bigWigToWig  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faAlign  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faCount  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faOneRecord  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSize  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/faSplit  ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/fetchChromSizes ./

RUN rsync -aP \
   rsync://hgdownload.soe.ucsc.edu/genome/admin/exe/linux.x86_64/wigToBigWig  ./

ENV PATH=${PATH}:/usr/src/kentUtils

# install SRA toolkit

WORKDIR $DEST

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

ENV PATH=${PATH}:/usr/src/picard/build/libs

# install biostar145820 (for shuffling bam files)
WORKDIR $DEST

RUN git clone --branch dev "https://github.com/lindenb/jvarkit.git"

WORKDIR $DEST/jvarkit

RUN $DEST/jvarkit/gradlew biostar145820

ENV PATH=${PATH}:/usr/src/jvarkit/dist

# reset working directory

WORKDIR /

# report done

CMD ["echo", "docker", "image", "built"]

WORKDIR /home
