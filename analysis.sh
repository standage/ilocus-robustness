#!/usr/bin/env bash
set -eo pipefail

download()
{
    # Download genome sequences/annotations, compute iLoci
    genhub-build.py --workdir=genomes \
                    --numprocs=4 \
                    --genome=Atha,Att6,Am32,Am10 \
                    download format prepare
}

prep()
{
    local speclabel=$1
    local oldlabel=$2
    local newlabel=$3

    mkdir -p $speclabel
    ln -fsn $(pwd)/genomes/${oldlabel}/${oldlabel}.iloci.gff3 \
            ${speclabel}/${newlabel}.iloci.gff3
    ln -fsn $(pwd)/genomes/${oldlabel}/${oldlabel}.iloci.fa \
            ${speclabel}/${newlabel}.iloci.fa
    ln -fsn $(pwd)/genomes/${oldlabel}/${oldlabel}.gdna.fa \
            ${speclabel}/${newlabel}.gdna.fa
}

mask_filter_repeats()
{
    local genome=$1
    local species=$2

    # Mask transposable elements and other repetitive DNA
    RepeatMasker -species $species -parallel $numprocs -frag 1000000 -lcambig \
                 -xsmall -gff ${genome}.gdna.fa \
      > ${genome}.rm.log 2>&1

    # Calculate the repetitive content of each iLocus
    bedtools coverage -a <(grep $'\tlocus\t' ${genome}.iloci.gff3) \
                      -b ${genome}.gdna.fa.out.gff \
        > ${genome}.iloci.coverage.tsv

    # Filter iLoci by repetitive content
    ./filter_iloci_by_coverage.py --perc 0.25 --bp 1000 \
                                  ${genome}.iloci.coverage.tsv \
        > ${genome}.iloci.filter-rm.txt \
        2> ${genome}.iloci.filter-rm.discard.txt
    ./select_seq.py ${genome}.iloci.filter-rm.txt ${genome}.iloci.fa \
        > ${genome}.iloci.filter-rm.fa
}


run_vmatch()
{
    local db=$1
    local query=$2
    local numprocs=$3

    # Create suffix tree
    mkvtree -db ${db}.gdna.fa.masked -indexname ${db}.gdna.fa.masked -dna \
            -pl 12 -allout -v

    # Split up iLocus input sequences
    gt splitfasta -numfiles 16 ${query}.iloci.filter-rm.fa

    # Search and report matches
    parallel --gnu --jobs $numprocs \
            vmatch -q ${query}.iloci.filter-rm.fa.{} -complete -e 5b \
                   -identity 95 -d -p -showdesc 0 ${db}.gdna.fa.masked \
                '>' ${db}.vmatch.{}.txt \
        ::: $(seq 1 $numprocs)

    # Compute mapping based on vmatch alignments
    ./ilocus_mapping.py --outfile=${db}.ilocus_map.txt \
                        --logfile=${db}.ilocus_map.log \
                        ${query}.iloci.gff3 \
                        ${db}.iloci.gff3 \
                        ${db}.vmatch.*.txt

    # Compute proportion of stable iLoci
    ilocus_count_txt=$(wc -l < ${query}.iloci.filter-rm.txt)
    ilocus_count_map=$(grep -v $'^None\t' ${db}.ilocus_map.txt | cut -f 1 | sort -u | wc -l)
    ilocus_count_stable=$(grep -v None ${db}.ilocus_map.txt | cut -f 1 | sort -u | wc -l)
    echo "iLocus count:"
    echo "    - total:    $ilocus_count_map"
    echo "    - filtered: $ilocus_count_txt"
    frac="$ilocus_count_stable / $ilocus_count_txt"
    echo "Stable iLoci: $frac = " $(Rscript -e "$frac" | cut -c 5-)
}


# Main method
download

prep Atha Att6 TAIR6
prep Atha Atha TAIR10
prep Amel Am10 OGS1.0
prep Amel Am32 OGS3.2

mask_filter_repeats Atha/TAIR6 16 &
mask_filter_repeats Amel/OGS1.0 16 &
wait

run_vmatch Atha/TAIR10 Atha/TAIR6 16 &
run_vmatch Amel/OGS3.2 Amel/OGS1.0 16 &
wait

