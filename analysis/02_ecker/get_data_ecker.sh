#!/bin/bash
##
## Retrieves GSE97179's (Ecker) data and sort of gets them into a tabular format
## with the tab-separated structure (no header):
##
## chr	pos	strand	mc_class	mc_count	total	methylated
##
## 14 Apr 2021
## Izaskun Mallona

NTHREADS=20
WD=/home/imallona/emanuel/GSE97179

mkdir -p "$WD"
cd "$WD"

# the annotation file for that accession
# got it from view-source:https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE97179

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97179/miniml/GSE97179_family.xml.tgz
tar xzvf *xml.tgz

raw="$(grep RAW GSE97179_family.xml)"

# ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE97nnn/GSE97179/suppl/GSE97179_RAW.tar
wget "$raw"

# that's a 683GB file
# let's get the individual files, no way of validating the md5sum (not available at GEO afaik)
tar -tf GSE97179_RAW.tar > children

# per cell/archived GEO unit CG-only cytosine report, gz compressed, tab separated
mkdir -p "$WD"/merged_cg

cat children | while read line
do
    echo "$line"
    tag=merged_"$(basename $line .tar.gz)".tsv.gz
    mkdir tmp_"$tag"; cd "$_"

    # that's to extract the cell/pool tarball
    tar xvf ../GSE97179_RAW.tar "$line"
    file "$line" # beware it's NOT a gz-compressed tar
    
    ## the strip-components is to flatten the directory structure
    tar xvf "$line" --strip-components 1  

    ## are there gz files?
    compressed=$(find . -name "*tsv.gz" | wc -l)
    ## are there nongz files?
    uncompressed=$(find . -name "*tsv" | wc -l)

    ## the GEO tarballs are weird and unconsistent, with sometimes compressed,
    # sometimes uncompressed cytosine reports
    if [ $compressed -eq 0 ] && [ $uncompressed != 0 ]
    then
        echo not compressed files
        grep -P -h "\tCG" *tsv | \
            pigz -p "$NTHREADS" > "$WD"/merged_cg/"$tag" 
    elif [ $uncompressed -eq 0 ] && [ $compressed != 0 ]
    then
        echo compressed files
        # uncompress, grep for CG context only, compress, send to file
        # there is no coordinate sorting, we rely on the archive's filenames
        pigz --decompress -p "$NTHREADS" --to-stdout *tsv.gz  | \
            grep -P "\tCG" | \
            pigz -p "$NTHREADS" > "$WD"/merged_cg/"$tag" 
    
    fi
     
    cd "$WD"
    if [ -d tmp_"$tag" ]; then rm -rf tmp_"$tag"; fi

done

rm children

# note file headers are:
# chr	pos	strand	mc_class	mc_count	total	methylated
