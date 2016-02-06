# Evaluating iLocus robustness across assembly versions

## Overview

Refinement of reference genome assemblies is an ongoing task for both mature model organisms and for novel research systems.
However, improved, more accurate assemblies come at the expense of disrupting the sequence-based coordinate system typically used for annotating the location of genome features.
Parsing an annotated genome into *iLoci* provides an alternative representation of the complete genome that is robust to assembly and annotation updates.

To evaluate this claim, we selected two model organisms (*Arabidopsis thaliana* and the honey bee *Apis mellifera*) for which multiple assemblies and annotations produced over the span of several years are available.
The 2005 TAIR6 release was the first annotation of the *A. thaliana* genome managed by The Arabidopsis Information Resource [TAIR][tair], while the 2010 TAIR10 release integrates TAIR's latest improvements to both the reference genome assembly and annotation [TAIR10][tair10].
For *A. mellifera*, the Honey Bee Genome Sequencing Consortium's assembly version 2 and the Official Gene Set 1 (OGS v1.0) constituted the initial description of the honey bee genome in 2006 [OGS1.0][ogs1.0], while assembly 4.5 and OGS v3.2 represent the consortium's latest improvements to the genome and corresponding annotation as of 2014 [OGS3.2][ogs3.2].

For these two species, we computed *iLoci* for both the earlier assembly/annotation (version *A*) and the newer assembly/annotation (version *B*) to evaluate the proportion of *iLoci* that remain stable between versions.
Here, we designate an *iLocus* as *stable* or if the following criteria are satisfied: the *iLocus* sequence from version *A* has a unique optimal global alignment to assembly version *B* with at least 95% identity; one *iLocus* from version *B* overlaps with at least 90% of the aligned *iLocus* sequence from version *A*; and the aligned *iLocus* sequence from *A* overlaps with at least 90% of that same *iLocus* from *B*.

[tair]: https://www.arabidopsis.org/
[tair10]: http://dx.doi.org/10.1093/nar/gkr1090
[ogs1.0]: http://dx.doi.org/10.1186/gb-2007-8-1-r13
[ogs3.2]: http://dx.doi.org/10.1186/1471-2164-15-86

## Procedure

Implementation is provided in `analysis.sh`.
Requires [GenHub][genhub] for downloading and preprocessing the data, and [vmatch][vmatch] for computing the alignments.
The recipes for TAIR6, OGS1.0, and OGS3.2 have not yet been integrated into GenHub's main development branch, so for now use the `robust` branch from GitHub.

```bash
time bash analysis.sh
```

[genhub]: https://github.com/standage/genhub
[vmatch]: http://www.vmatch.de
