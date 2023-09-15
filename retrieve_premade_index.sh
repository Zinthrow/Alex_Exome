wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz
mv GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz data/
tar -xzf data/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index.tar.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
gunzip GCA_000001405.15_GRCh38_full_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_full_analysis_set.fna data/GCA_000001405.15_GRCh38_full_analysis_set.fna
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_full_analysis_set.fna.fai
mv GCA_000001405.15_GRCh38_full_analysis_set.fna.fai data/GCA_000001405.15_GRCh38_full_analysis_set.fna.fai