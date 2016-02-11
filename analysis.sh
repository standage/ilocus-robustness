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
    # Create clean working directories for this analysis
    mkdir -p Atha/
    mkdir -p Amel/

    # iLocus sequences (older assembly/annotation version)
    ln -fsn $(pwd)/genomes/Att6/Att6.iloci.fa Atha/TAIR6.iloci.fa
    ln -fsn $(pwd)/genomes/Am10/Am10.iloci.fa Amel/OGS1.0.iloci.fa

    # Whole genome sequences (newer assembly version)
    ln -fsn $(pwd)/genomes/Atha/Atha.gdna.fa Atha/TAIR10.gdna.fa
    ln -fsn $(pwd)/genomes/Am32/Am32.gdna.fa Amel/OGS3.2.gdna.fa

    # iLocus annotations (older version)
    ln -fsn $(pwd)/genomes/Att6/Att6.iloci.gff3 Atha/TAIR6.iloci.gff3
    ln -fsn $(pwd)/genomes/Am10/Am10.iloci.gff3 Amel/OGS1.0.iloci.gff3

    # iLocus annotations (newer version)
    ln -fsn $(pwd)/genomes/Atha/Atha.iloci.gff3 Atha/TAIR10.iloci.gff3
    ln -fsn $(pwd)/genomes/Am32/Am32.iloci.gff3 Amel/OGS3.2.iloci.gff3
}

run_vmatch()
{
    local db=$1
    local query=$2

    # Create suffix tree
    mkvtree -db ${db}.gdna.fa -indexname ${db}.gdna.fa -dna -pl 12 -allout -v

    # Search and report matches
    vmatch -q ${query}.iloci.fa -mum -l 400 -d -p -identity 50 -showdesc 0 \
           ${db}.gdna.fa \
        > ${db}.vmatch.txt
    ./ilocus_mapping.py ${query}.iloci.gff3 \
                        ${db}.iloci.gff3 \
                        ${db}.vmatch.txt \
        > ${db}.ilocus_map.txt
}

download
prep
run_vmatch Atha/TAIR10 Atha/TAIR6 &
run_vmatch Amel/OGS3.2 Amel/OGS1.0 &
wait

