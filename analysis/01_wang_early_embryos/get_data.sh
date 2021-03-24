#!/bin/bash
##
## Downloads data from https://www.nature.com/articles/s41467-021-21409-8#data-availability
##  Wang et al CpG/GpC meth from GEO GSE136718 and converts it to bedgraph
##  and prettifies them for downstream analysis
##
## Example output for a given cell (CpG coordinate, beta-value, tab-separated):
##  chr1:3000121    1
##  chr1:3000180    0.25
##  chr1:3000666    0
##  chr1:3000699    1
##
## Izaskun Mallona
## 24 March 2021

NTHREADS=32
bw2bed=~/soft/kent/userApps.v390/bin/bigWigToBedGraph
GSE=GSE136718

cd ~/smp_test/

wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE136nnn/GSE136718/suppl//GSE136718_RAW.tar" \
     --output-file "$GSE".tar

tar xvf "$GSE".tar

# get plain text files from bigwigs, with just two columns:
# chr:start \t beta_value
# reason why:
#   we need the chr:start as identifier, the end coordinate we don't care
#   it's CG methylation, so it's start +1, and the beta value is the variable of interest
# processed in blocks of $NTHREADS files (input bigwigs) each time

N=$NTHREADS

(
    for fn in $(find . -name "*WCG*bw"| xargs -n"$NTHREADS")
    do 
        ((i=i%N)); ((i++==0)) && wait
        curr=$(basename "$fn" .bw)
        echo $curr $i
        
        "$bw2bed" "$fn" stdout | awk '{OFS=FS="\t"; print $1":"$2, $4}' > "$curr".dat &
    done
)

# how many CpGs each?

find . -name "*WCG*dat" -type f -exec wc -l {} \; | sort -k1 -V -r

# compress & clean

pigz -p "$NTHREADS" *dat
rm ./*bw
rm ./"$GSE".tar
