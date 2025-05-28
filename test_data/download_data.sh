#!/bin/sh
# Download data for genotype-based deconvolution methods (popscle tutorial dataset)
outputdir="test_data"
mkdir -p $outputdir && cd $outputdir
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1leRqeX9BDvZHtSncG6QTbQmUsQJOPlZq' -O jurkat_293t_exons_only.vcf.withAF.vcf.gz
wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=12iyB7J55wRtC01L7dYuPlQKzcQdjZKe6' -O jurkat_293t_downsampled_n500_full_bam.bam.bai

wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=14HZgDCn1wBjDMDrQNIWiegA2rgu-Efz5' -O jurkat_293t_demuxlet.best

# need to download this with gdown because it's a big file and it requires a confirmation dialogue
FILEID="16I0P2eWj9tqkJ9BwakcZ2RMN_VCXkAv1"
FILENAME="jurkat_293t_downsampled_n500_full_bam.bam"
pip install gdown
gdown --id $FILEID -O $FILENAME

# block this part
#wget --no-check-certificate 'https://docs.google.com/uc?export=download&id=1MmwEiOsdzEfRdXS6oXXBwMJXUovKWcni' -O final_res.zip
#unzip final_res.zip
#rm final_res.zip
#mv final_res/jurkat_293t_demuxlet.best .
#rm -rf final_res

# To run souporcell, unzip VCF file
gzip -dk jurkat_293t_exons_only.vcf.withAF.vcf.gz

# Download subset reference genome
wget --no-check-certificate https://figshare.com/ndownloader/files/43102459 -O genome_chr1.fa
wget --no-check-certificate https://figshare.com/ndownloader/files/43102453 -O genome_chr1.fa.fai
# source: http://cf.10xgenomics.com/supp/cell-exp/refdata-cellranger-hg19-3.0.0.tar.gz

# Download common variants
wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1lw4T6d7uXsm9dt39ZtEwpuB2VTY3wK1y" -O common_variants_hg19.vcf && rm -rf /tmp/cookies.txt
wget https://master.dl.sourceforge.net/project/cellsnp/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg19.vcf.gz
# Processed by bcftools query -f '%CHROM:%POS\n' common_variants_hg19.vcf > common_variants_hg19_list.vcf
wget https://s3-eu-west-1.amazonaws.com/pfigshare-u-files/41773431/common_variants_hg19_list.vcf?X-Amz-Algorithm=AWS4-HMAC-SHA256 &
X-Amz-Credential=AKIAIYCQYOYV5JSSROOA/20230731/eu-west-1/s3/aws4_request &
X-Amz-Date=20230731T153655Z &
X-Amz-Expires=10 &
X-Amz-SignedHeaders=host &
X-Amz-Signature=283575c6ccb3104b8b95684e6d955abd28b47db71c18d1eeec99ae5dab65ff7b

# Download simualated data
# Can also run Rscript simulation.R
# Barcodes
wget --no-check-certificate https://figshare.com/ndownloader/files/41773428 -O barcodes.tsv
# common variant list for scSplit
wget --no-check-certificate https://figshare.com/ndownloader/files/41773431 -O common_variants_hg19_list.vcf
# simulated HTO and RNA counts
wget --no-check-certificate https://figshare.com/ndownloader/files/41779407 -O hto.zip
unzip hto.zip
wget --no-check-certificate https://figshare.com/ndownloader/files/41779260 -O rna.zip
unzip rna.zip
rm hto.zip
rm rna.zip
