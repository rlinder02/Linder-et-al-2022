FROM bioconductor/bioconductor_docker:RELEASE_3_15

COPY conda.yml .

RUN R -e 'BiocManager::install("BSgenome.Scerevisiae.UCSC.sacCer3")'
RUN R -e 'BiocManager::install("DECIPHER")'
RUN R -e 'BiocManager::install("DEGreport")'
RUN R -e 'BiocManager::install("EnsDb.Scerevisiae.v92")'
RUN R -e 'BiocManager::install("GO.db")'
RUN R -e 'BiocManager::install("Gviz")'
RUN R -e 'BiocManager::install("KEGG.db")'
RUN R -e 'BiocManager::install("KEGGREST")'
RUN R -e 'BiocManager::install("motifbreakR")'
RUN R -e 'BiocManager::install("MotifDb")'
RUN R -e 'BiocManager::install("org.Sc.sgd.db")'
RUN R -e 'BiocManager::install("pathview")'
RUN R -e 'BiocManager::install("phangorn")'
RUN R -e 'BiocManager::install("quadprog")'
RUN R -e 'BiocManager::install("RIPSeeker")'
RUN R -e 'BiocManager::install("seqinr")'
RUN R -e 'BiocManager::install("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")'

RUN R -e "install.packages('abind',dependencies=TRUE)"
RUN R -e "install.packages('colorspace',dependencies=TRUE)" 
RUN R -e "install.packages('cowplot',dependencies=TRUE)" 
RUN R -e "install.packages('data.table',dependencies=TRUE)" 
RUN R -e "install.packages('DescTools',dependencies=TRUE)"
RUN R -e "install.packages('devtools',dependencies=TRUE)"
RUN R -e "install.packages('doParallel',dependencies=TRUE)" 
RUN R -e "install.packages('egg',dependencies=TRUE)"
RUN R -e "install.packages('formattable',dependencies=TRUE)" 
RUN R -e "install.packages('GGally',dependencies=TRUE)"
RUN R -e "install.packages('ggcorrplot',dependencies=TRUE)"
RUN R -e "install.packages('ggdendro',dependencies=TRUE)"
RUN R -e "install.packages('ggExtra',dependencies=TRUE)"
RUN R -e "install.packages('ggh4x',dependencies=TRUE)"
RUN R -e "install.packages('ggpbeeswarm',dependencies=TRUE)"
RUN R -e "install.packages('ggpmisc',dependencies=TRUE)"
RUN R -e "install.packages('ggpubr',dependencies=TRUE)"
RUN R -e "install.packages('GISTools',dependencies=TRUE)"
RUN R -e "install.packages('gprofiler2',dependencies=TRUE)"
RUN R -e "install.packages('gridGaphics',dependencies=TRUE)"
RUN R -e "install.packages('here',dependencies=TRUE)" 
RUN R -e "install.packages('hexbin',dependencies=TRUE)"
RUN R -e "install.packages('HMM',dependencies=TRUE)"
RUN R -e "install.packages('lemon',dependencies=TRUE)"
RUN R -e "install.packages('magick',dependencies=TRUE)"
RUN R -e "install.packages('NbClust',dependencies=TRUE)"
RUN R -e "install.packages('pcaPP',dependencies=TRUE)"
RUN R -e "install.packages('plotrix',dependencies=TRUE)"
RUN R -e "install.packages('qqman',dependencies=TRUE)"
RUN R -e "install.packages('RSQLite',dependencies=TRUE)"
RUN R -e "install.packages('R.utils',dependencies=TRUE)" 
RUN R -e "install.packages('svMisc',dependencies=TRUE)"  
RUN R -e "install.packages('tictoc',dependencies=TRUE)"
RUN R -e "install.packages('tidyverse',dependencies=TRUE)" 
RUN R -e "install.packages('vcfR',dependencies=TRUE)" 

RUN R -e "install.packages('rBLAST', repos = 'https://mhahsler.r-universe.dev')"

RUN curl -sSL https://repo.continuum.io/miniconda/Miniconda3-py38_4.12.0-Linux-x86_64.sh -o /tmp/miniconda.sh \
    && bash /tmp/miniconda.sh -bfp /usr/local/ \
    && rm -rf /tmp/miniconda.sh

RUN conda env update --name root -f conda.yml && conda clean -a
