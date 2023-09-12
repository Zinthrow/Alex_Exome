FROM ubuntu:20.04

ENV DEBIAN_FRONTEND=noninteractive

# Install basic tools
RUN apt update && apt install -y wget unzip default-jdk bwa samtools bowtie2 fastqc

WORKDIR /tools
# FASTQC and Trimmomatic quality control tools
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip && \
    cp Trimmomatic-0.39/trimmomatic-0.39.jar /tools/ && \
    echo 'alias trimmomatic="java -jar /tools/trimmomatic-0.39.jar"' >> ~/.bashrc && \
    chmod 777 /tools/trimmomatic-0.39.jar && \
    rm -r Trimmomatic-0.39 Trimmomatic-0.39.zip

# Install GATK, or just use gatks offically created image
#RUN wget -O gatk.tar.bz2 https://github.com/broadinstitute/gatk/releases/download/4.1.9.0/gatk-4.1.9.0.tar.bz2 && \
#    tar xjf gatk.tar.bz2 && \
#    rm gatk.tar.bz2

# propriatary software must be obtained from qiagen
COPY  annovar.latest.tar.gz .
RUN tar xzf annovar.latest.tar.gz && \
     rm annovar.latest.tar.gz

# Add tools to PATH
#ENV PATH="/tools/gatk-4.1.9.0/:$PATH"
ENV PATH="/tools/annovar/:$PATH" 

CMD ["/bin/bash"]
