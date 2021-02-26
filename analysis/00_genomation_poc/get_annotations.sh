#!/bin/bash
##
## Fetches several annotations from UCSC, merging the overlapping/bookended ones
##   and reporting them as bed files
##
## Caution 0/1 counting coordinates not harmonized!
##
## 25 February 2021
## Izaskun Mallona

# pending cpgislands, pmds, other entropy scores

set -o errexit
set -o pipefail

BEDTOOLS=~/soft/bedtools/bedtools-2.29.2/bin/bedtools
KENT=~/soft/kent/userApps.v390/bin/
## ~/soft/bedtools/bedtools2/bin/bedtools --version
## bedtools v2.27.1

WD="$HOME/smp_test"
NTHREADS=10

mkdir -p "$WD"/data
cd $_

## base annotation start

mysql --user=genome --host=genome-mysql.cse.ucsc.edu -A -N \
      -e "select chrom, size from mm10.chromInfo" > mm10.genome

wget http://hgdownload.soe.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz

pigz -d -p "$NTHREADS" mm9ToMm10.over.chain.gz

## base annotation end

## promoters start

## THIS IS LIKELY WRONG
mkdir -p $WD/data/promoters
cd $_

mysql --user=genome \
      --host=genome-mysql.cse.ucsc.edu -A \
      -B -N \
      -e "select 
           chrom, 
           txStart,
           txEnd,
           name,
           strand           
          from 
           mm10.refGene;" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4,"0",$5}' | \
    $BEDTOOLS sort > refGene_mm10.bed

## bedtools expanding overlapping transcripts and getting the upstream region,
## taking into account the strand
## this is most likely wrong

$BEDTOOLS merge -s -i refGene_mm10.bed -c 4,6 -o distinct | \
    $BEDTOOLS slop -i - -g ../mm10.genome -l 1000 -r 500 -s > mm10_nonoverlapping_promoters.bed

pigz -p $NTHREADS *.bed

## promoters end

## lamina start

mkdir $WD/data/lamina
cd $_

mysql --user=genome --host=genome-euro-mysql.soe.ucsc.edu -A \
      -B -N -e \
      "select chrom,
              chromStart,
              chromEnd,
              sumData
       from mm9.laminB1_AC" | \
    awk '{OFS=FS="\t"; print $1,$2,$3,$4}' | \
    $BEDTOOLS sort > mm9_laminB1_AC.bed

"$KENT"/liftOver mm9_laminB1_AC.bed \
       ../mm9ToMm10.over.chain \
       mm9_laminB1_AC_liftovered_mm10.bed \
       mm9_laminB1_AC_liftovered_mm10_unmapped.bed


## maybe weird, but what about collapsing regions with positive
##  values if they're closer than 5kbp?
$BEDTOOLS sort -i mm9_laminB1_AC_liftovered_mm10.bed |\
    awk 'BEGIN {FS=OFS="\t" }; { 
      if($4 < 0) print $1,$2,$3,"n", "0", "+"; 
      else print $1, $2,$3,"p","0","-"
      }'  | \
        $BEDTOOLS merge -d 5000 -c 4 -o distinct -s > lamina_mm10_test.bed

pigz -p $NTHREADS *.bed

## lamina end

## replichip start

mkdir $WD/data/replichip
cd $_

# data from
# https://www.encodeproject.org/experiments/ENCSR000AXO/

## replicate 1
wget https://www.encodeproject.org/files/ENCFF001JUP/@@download/ENCFF001JUP.bigWig

# replicate 2
wget https://www.encodeproject.org/files/ENCFF001JUQ/@@download/ENCFF001JUQ.bigWig

for bw in $(find . -name "*bigWig")
do
    
    "$KENT"/bigWigToBedGraph "$bw" "$(basename $bw .bigWig)".bed
    
    "$KENT"/liftOver "$(basename $bw .bigWig)".bed \
           ../mm9ToMm10.over.chain \
           "$(basename $bw .bigWig)"_liftovered_mm10.bed \
           "$(basename $bw .bigWig)"_liftovered_mm10_unmapped.bed
done

pigz -p "$NTHREADS" *bed

## replichip end

## chromHMM start

## chromHMM end

## windows start
mkdir -p $WD/data/windows
cd $_

$BEDTOOLS makewindows -g ../mm10.genome -w 100000 > mm10_100k_bp_windows.bed
$BEDTOOLS makewindows -g ../mm10.genome -w 10000 > mm10_10k_bp_windows.bed
$BEDTOOLS makewindows -g ../mm10.genome -w 1000 > mm10_1k_bp_windows.bed

pigz -p "$NTHREADS" *bed

## windows end

## cpgislands start

## notice we keep islands overlapping repeats
mkdir -p "$WD"/data/cpgIslandExtUnmasked

cd $_

mysql --user=genome --host=genome-euro-mysql.soe.ucsc.edu -A \
      -B -N -e \
      "select chrom,
              chromStart,
              chromEnd,
              obsExp
       from mm10.cpgIslandExtUnmasked" | \
    $BEDTOOLS sort > mm10_cpgIslandExtUnmasked.bed

pigz -p "$NTHREADS" *bed

## cpgislands end


## chromHMM ESC mm10 start

mkdir -p "$WD"/data/chromHMM
cd $_

wget https://github.com/guifengwei/ChromHMM_mESC_mm10/raw/master/mESC_E14_12_dense.annotated.bed.gz

zcat mESC_E14_12_dense.annotated.bed.gz | sed '1d' | cut -f1-4 > mESC_E14_12_chromHMM.bed

rm mESC_E14_12_dense.annotated.bed.gz
pigz -p "$NTHREADS" *bed

## chromHMM ESC mm10 end

# compress on exit

cd $WD/data

pigz -p "$NTHREADS" mm9ToMm10.over.chain

# maybe delete bigWigs etc
