# All the UCSC scripts can be found here: https://github.com/ENCODE-DCC/kentUtils
UCSC_scripts=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/tools/kentUtils/bin/linux.x86_64/

# The liftover files found here: https://hgdownload.soe.ucsc.edu/gbdb/hg19/liftOver/
liftover_file=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/liftover_ucsc/liftOver/hg19ToHg38.over.chain

# The repliseq bigwig files are here: https://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeFsuRepliChip/
repliseq_dir=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/repliseq_data/wgEncodeFsuRepliChip

# Getting the chrom sizes at: https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/
chrom_sizes=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/liftover_ucsc/hg38.chrom.sizes

output_dir=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/repliseq_data/hg38_repliseq_liftover
tmp_dir=/home/anna_h/bioinf_isilon/Research/CRESSWELL/Internal/ahakobyan/tmp

# mkdir $output_dir

for bigwig_file in `ls $repliseq_dir/*.bigWig`; do
    echo ${bigwig_file}
    bigwig_filename=$(basename ${bigwig_file})
    bedgraph_filename=${bigwig_filename%.bigWig}".bedGraph"
    bedgraph_liftover=${bigwig_filename%.bigWig}"_hg38.bedGraph"

    bigwig_hg38_filename=${bigwig_filename%.bigWig}"_hg38.bigWig"

    # ${UCSC_scripts}/bigWigToBedGraph ${bigwig_file} ${tmp_dir}/${bedgraph_filename}

    # ${UCSC_scripts}/liftOver <(sort -k1,1 -k2,2n ${tmp_dir}/${bedgraph_filename}) ${liftover_file} ${tmp_dir}/${bedgraph_liftover} ${tmp_dir}/${bedgraph_liftover}".unMapped"

    ### the awk code skips all overlapping regions and prints only the first one
    sort -k1,1 -k2,2n ${tmp_dir}/${bedgraph_liftover} | \
            awk 'BEGIN{OFS="\t"}
            NR == 1 {
                chr=$1; # reading the values
                start=$2;
                end=$3;
                val=$4;
                next;
            }
            {
                if (chr != $1 || end <= $2) {
                    print chr, start, end, val;
                    chr = $1;
                    start = $2;
                    end = $3;
                    val = $4;
                } else {
                    next
                }
            } END {
                print chr, start, end, val;
            }' > ${tmp_dir}/${bedgraph_liftover}"._sorted.bedGraph"

    
    $UCSC_scripts/bedGraphToBigWig ${tmp_dir}/${bedgraph_liftover}"._sorted.bedGraph" ${chrom_sizes} ${output_dir}/${bigwig_hg38_filename}

done


