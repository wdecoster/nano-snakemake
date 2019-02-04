# nano-snakemake

WORK IN PROGRESS

A set of [snakemake](https://snakemake.readthedocs.io/en/stable/) pipelines for nanopore whole genome sequencing processing and analysis.


- SV analysis pipeline using:
  - [ngmlr](https://github.com/philres/ngmlr)
  - [sniffles](https://github.com/fritzsedlazeck/Sniffles)
  - [NanoSV](https://github.com/mroosmalen/nanosv)
  - [survivor](https://github.com/fritzsedlazeck/SURVIVOR)
  - [mosdepth](https://github.com/brentp/mosdepth)
  - [minimap2](https://github.com/lh3/minimap2)
  - [cyvcf2](https://github.com/brentp/cyvcf2)
  - [vcfanno](https://github.com/brentp/vcfanno)
  - [vcftools](https://vcftools.github.io/index.html)
  - [samtools](https://github.com/samtools/samtools)
  - [matplotlib](https://github.com/matplotlib/matplotlib)
  - [seaborn](https://github.com/mwaskom/seaborn)

These dependencies can be installed using `conda create -f environment.yaml`


The folder "scripts" contains scripts necessary for the pipeline.   
The folder "extra_scripts" contains scripts I wrote while playing around with the results, mainly visualizations.
