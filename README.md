# Evaluating iLocus robustness across assembly versions

In this analysis I am evaluating the robustness of *iLoci* between two assembly/annotation versions of the same genome.
The precise definition of an *iLocus*, described in detail elsewhere, is outside the scope of this analysis.
Suffice it to say that an *iLocus* is a small region of the genome encapsulating a single gene, a set of overlapping genes, or an intergenic region.
It's important to note that *iLoci* typically do not overlap, although in some cases overlap is allowed.

## Input data

The input data for this analysis consists of genome assemblies (sequences in Fasta format) and associated annotations (genome features in GFF3 format).
We are looking at data from two species, the plant *Arabidopsis thaliana* and the honey bee *Apis mellifera*.
For each species, we are comparing data from one assembly/annotation version to another.

- *Arabidopsis thaliana* (Atha)
    - TAIR6
    - TAIR10
- *Apis mellifera* (Amel)
    - Amel 2.0 / OGS 1.0
    - Amel 4.5 / OGS 3.2

Each of these four assembly/annotation pairs is processed to compute iLoci.
The details of this processing is also outside the scope of this analysis, but commands for doing the processing are provided in the **Obtaining input data** section.

### Obtaining input data

The genome assemblies and annotations needed for this analysis are available from public databases.
We use GenHub to download the data files, pre-process them into a consistent format, compute iLoci, and extract iLocus sequences.
Please see [PREREQS.md](PREREQS.md) for technical details on installing GenHub and its dependencies.

```bash
genhub-build.py --workdir=genomes --numprocs=4 \
                --genome=Atha,Att6,Am32,Am10 \
                download format prepare
```

### Setting up a clean working directory

Assuming everything ran correctly, GenHub created the `genomes/` directory and subdirectories for each assembly/annotation version.
These subdirectories contain quite a few files, including ancillary files that are irrelevant to this analysis.
Let's make a clean directory---to isolate the work done in this analysis from the processing done by GenHub---and create links to the files we actually need.
For each annotated genome we need:
- iLocus annotation
- iLocus sequences
- whole genome sequences

```bash
mkdir Atha/
mkdir Amel/

# TAIR6
ln -s $(pwd)/genomes/Att6/Att6.iloci.gff3 Atha/TAIR6.iloci.gff3
ln -s $(pwd)/genomes/Att6/Att6.iloci.fa Atha/TAIR6.iloci.fa
ln -s $(pwd)/genomes/Att6/Att6.gdna.fa Atha/TAIR6.gdna.fa

# TAIR10
ln -s $(pwd)/genomes/Atha/Atha.iloci.gff3 Atha/TAIR10.iloci.gff3
ln -s $(pwd)/genomes/Atha/Atha.iloci.fa Atha/TAIR10.iloci.fa
ln -s $(pwd)/genomes/Atha/Atha.gdna.fa Atha/TAIR10.gdna.fa

# Amel 2.0 / OGS 1.0
ln -s $(pwd)/genomes/Am10/Am10.iloci.gff3 Amel/OGS1.0.iloci.gff3
ln -s $(pwd)/genomes/Am10/Am10.iloci.fa Amel/OGS1.0.iloci.fa
ln -s $(pwd)/genomes/Am10/Am10.gdna.fa Amel/OGS1.0.gdna.fa

# Amel 4.5 / OGS 3.2
ln -s $(pwd)/genomes/Am32/Am32.iloci.gff3 Amel/OGS3.2.iloci.gff3
ln -s $(pwd)/genomes/Am32/Am32.iloci.fa Amel/OGS3.2.iloci.fa
ln -s $(pwd)/genomes/Am32/Am32.gdna.fa Amel/OGS3.2.gdna.fa
```

## iLocus alignment

So far the goal of this analysis has only been described in vague terms: we want to determine how robust or stable iLoci are from one assembly/annotation version to another.
Now we need to be more precise about what we mean by *stable*.
An iLocus *i<sub>A</sub>* from assembly/annotation version *A* is __*stable*__ if it satisfies the following criteria.
- at least 90% of *i<sub>A</sub>* can be aligned to assembly *B*
- the aligned (sub)sequences have at least 95% identity
- there exists at least one locus *i<sub>B</sub>* from asssembly/annotation version *B* that has at least 90% reciprocal overlap with *i<sub>A</sub>*; that is, *i<sub>A</sub>* overlaps with at least 90% of *i<sub>B</sub>*, and *i<sub>B</sub>* overlaps with at least 90% of *i<sub>A</sub>*

We can use the vmatch software suite to compute the alignments.
However, since we are computing approximate matches instead of exact matches, highly repetitive sequences have the potential to drastically increase the runtime required to compute the alignments.
So we will want to annotate repetitive sequences on genome *A* and filter out corresponding iLoci with repetitive content above some threshold.
The iLoci that satisfy the filtering criteria will be aligned to genome *B*, and the proportion of these iLoci that align (as opposed to all iLoci) will be used to assess iLocus stability.

### Annotation of repeats

We will use RepeatMasker to annotate transposable elements and other repetitive DNA in the genome.
We will then use bedtools to compute the percentage of each iLocus occupied by repeats, and custom scripts to filter out iLoci with repeats occupying at least 25% of their length or 500 bp total.

```bash
#----- Arabidopsis thaliana -----
RepeatMasker -species viridiplantae -parallel 16 -frag 1000000 -lcambig \
             -xsmall -gff Atha/TAIR6.gdna.fa \
    > Atha/TAIR6.rm.log 2>&1

bedtools coverage -a <(grep $'\tlocus\t' Atha/TAIR6.iloci.gff3) \
                  -b Atha/TAIR6.gdna.fa.out.gff \
    > Atha/TAIR6.iloci.coverage.tsv

./ilocus_filter.py --perc 0.25 --bp 500 Atha/TAIR6.iloci.coverage.tsv \
    > Atha/TAIR6.iloci.filter-rm.txt \
    2> Atha/TAIR6.iloci.filter-rm.discard.txt

./select_seq.py Atha/TAIR6.iloci.filter-rm.txt Atha/TAIR6.iloci.fa \
    > Atha/TAIR6.iloci.filter-rm.fa

#----- Apis mellifera -----
RepeatMasker -species insects -parallel 16 -frag 1000000 -lcambig \
             -xsmall -gff Amel/OGS1.0.gdna.fa \
    > Amel/OGS1.0.rm.log 2>&1

bedtools coverage -a <(grep $'\tlocus\t' Amel/OGS1.0.iloci.gff3) \
                  -b Amel/OGS1.0.gdna.fa.out.gff \
    > Amel/OGS1.0.iloci.coverage.tsv

./ilocus_filter.py --perc 0.25 --bp 500 Amel/OGS1.0.iloci.coverage.tsv \
    > Amel/OGS1.0.iloci.filter-rm.txt \
    2> Amel/OGS1.0.iloci.filter-rm.discard.txt

./select_seq.py Amel/OGS1.0.iloci.filter-rm.txt Amel/OGS1.0.iloci.fa \
    > Amel/OGS1.0.iloci.filter-rm.fa
```

### Computing alignments

Having selected iLoci that match our filtering criteria, we are ready to compute the alignments.
First we must create an index of the genome sequences to facilitate rapid searches.

```bash
mkvtree -db Atha/TAIR10.gdna.fa -indexname Atha/TAIR10.gdna.fa -dna -pl 12 -allout -v
mkvtree -db Am32/OGS3.2.gdna.fa -indexname Atha/OGS3.2.gdna.fa -dna -pl 12 -allout -v
```

And now we can compute the alignments.
With the `vmatch` command, we can only enforce the `95%` identity criterion.
The other criteria regarding receiprocal overlap will have to be handled with post-processing.

```bash
vmatch -q Atha/TAIR6.iloci.filter-rm.fa -complete -e 5b -identity 95 -d -p -showdesc 0 Atha/TAIR10.gdna.fa \
    > Atha/TAIR10.vmatch.txt
vmatch -q Amel/OGS1.0.iloci.filter-rm.fa -complete -e 5b -identity 95 -d -p -showdesc 0 Amel/OGS3.2.gdna.fa \
    > Amel/OGS3.2.vmatch.txt
```

### Filtering alignments

Now that we have filtered iLoci based on repeat content, aligned those iLoci from assembly *A* to assembly *B*, we need to check how these alignments line up with iLoci in the new assembly.
The `ilocus_mapping.py` script reads the vmatch output, along with the iLocus annotations from both assemblies, and determines which iLoci are mapped to the new assembly and iLocus annotation according to the overlap criteria specified above.

```bash
./ilocus_mapping.py --outfile=Atha/TAIR10.ilocus_map.txt --logfile=Atha/TAIR10.ilocus_map.log \
    Atha/TAIR6.iloci.gff3 Atha/TAIR10.iloci.gff3 Atha/TAIR10.vmatch.txt
./ilocus_mapping.py --outfile=Amel/OGS3.2.ilocus_map.txt --logfile=Amel/OGS3.2.ilocus_map.log \
    Amel/OGS1.0.iloci.gff3 Amel/OGS3.2.iloci.gff3 Atha/OGS3.2.vmatch.txt
```

## Results

### *Arabidopsis thaliana*

```
===== iLocus stability report (A=Atha/TAIR6, B=Atha/TAIR10) =====
    - total iLoci from A: 37299
    - filtered iLoci:     33730
    - aligned iLoci:      33717 / 33730 =  0.9996146
    - stable iLoci:       28890 / 33730 =  0.8565076

===== Breakdown: stable iLoci =====
fiLocus: 4
niLocus: 730
piLocus: 22749
iiLocus: 5395
ciLocus: 12

===== Breakdown: non-stable iLoci =====
niLocus: 39
fiLocus: 5
piLocus: 3758
iiLocus: 4602
ciLocus: 5
```

### *Apis mellifera*

```
===== iLocus stability report (A=Amel/OGS1.0, B=Amel/OGS3.2) =====
    - total iLoci from A: 22784
    - filtered iLoci:     19447
    - aligned iLoci:      18010 / 19447 =  0.9261069
    - stable iLoci:       5519 / 19447 =  0.283797

===== Breakdown: stable iLoci =====
fiLocus: 2247
piLocus: 2697
iiLocus: 567
ciLocus: 8

===== Breakdown: non-stable iLoci =====
fiLocus: 5743
piLocus: 7201
iiLocus: 4223
ciLocus: 98
```

## Appendix: motivation

Refinement of reference genome assemblies is an ongoing task for both mature model organisms and for novel research systems.
Annotated genome features are almost always described by their absolute position on an assembled genomic sequence.
As a result, if subsequent refinements to an assembly add a nucleotide here and remove a structural variant there, the whole coordinate system upon which the annotation relies is changed.
This is why gene features annotated on one version of a genome assembly (for example, hg18) are invalid on another version (hg19).

However, there will be a substantial portion of the genome for which the local genomic neighborhood remains unchanged between assembly versions.
Maybe a particular gene starts at position 568,366 of the scaffold, whereas before it started at position 568,270.
But if the sequence of the gene itself (along with a small amount of flanking intergenic space) is the same, that gene annotation is stable, even though we have to specify its location using different coordinates in the two assemblies.

This is one of the primary motivations for parsing an annotated genome sequence into iLoci, and using these iLoci as the coordinate system in favor of entire chromosome or scaffold sequences.
Genes or intergenic regions whose genomic context remains unchanged will be identical between assembly/annotation versions, and while their global coordinates may change their coordinates relative to the iLocus will be stable.
The extent to which iLoci are stable across assembly/annotation versions is what this analysis seeks to quantify.
