#!/usr/bin/env bash
set -eo pipefail

download()
{
    # Download genome sequences/annotations, compute iLoci
    genhub-build.py --workdir=genomes \
                    --numprocs=4 \
                    --genome=Atha,Att6,Am32,Am10 \
                    download format prepare stats
}

prep()
{
    # Create clean working directories for this analysis
    mkdir -p Atha
    mkdir -p Amel
    ln -fsn $(pwd)/genomes/Atha/Atha.iloci.fa Atha/TAIR10.iloci.fa
    ln -fsn $(pwd)/genomes/Att6/Att6.iloci.fa Atha/TAIR6.iloci.fa
    ln -fsn $(pwd)/genomes/Am32/Am32.iloci.fa Amel/OGS3.2.iloci.fa
    ln -fsn $(pwd)/genomes/Am10/Am10.iloci.fa Amel/OGS1.0.iloci.fa
}

run_vmatch()
{
    local db=$1
    local query=$2
    # Create suffix tree
    mkvtree -db $db -indexname $db -dna -pl 12 -allout -v \
        > ${db}.mkvtree.log 2>&1
    # Search and report matches
    vmatch -q $query -complete -d -p -identity 95 -showdesc 0 $db \
        > ${db}.vmatch.out.txt \
        2> ${db}.vmatch.log
}

download
prep
run_vmatch Atha/TAIR10.iloci.fa Atha/TAIR6.iloci.fa
run_vmatch Amel/OGS3.2.iloci.fa Amel/OGS1.0.iloci.fa
