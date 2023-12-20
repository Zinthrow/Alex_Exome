reads=`pwd`"/data"
base_name="10847101"

# download clinvar make database using
# git clone https://github.com/mobidic/update_annovar_db.git

# you'll need to run this using your own paths and naming conventions
# sudo python3 update_resources.py -d clinvar -hp /home/alarsen/Projects/Alex_Exome/data/humandb -a /home/alarsen/Projects/Alex_Exome/data/2023Dec -g GRCh38
# I had trouble with this / only attempted it once since this package seems to prefer a non-containerized version of annovar which I didn't see the point finagling

# using the ftp site I downloaded clinvar instead
# wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20231217.vcf.gz
# wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20231217.vcf.gz.tbi

docker run -v $reads:$reads annovar_bioinformatics bash -c "
convert2annovar.pl -format vcf4 clinvar_20231217.vcf.gz > humandb/hg38_clinvar20231217.txt"

# run clinvar
docker run -v $reads:$reads annovar_bioinformatics bash -c "
table_annovar.pl $reads/${base_name}_variants.vcf $reads/humandb/ -buildver hg38 -out $reads/${base_name}_clinvar -remove -protocol clinvar20231217 -operation g,f -nastring . -vcfinput"
