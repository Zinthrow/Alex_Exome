wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
gunzip GCF_000001405.40_GRCh38.p14_genomic.fna.gz
mv GCF_000001405.40_GRCh38.p14_genomic.fna data/
# bwa index -a bwtsw data/GCF_000001405.40_GRCh38.p14_genomic.fna
bowtie2-build data/GCF_000001405.40_GRCh38.p14_genomic.fna GRCh38_index