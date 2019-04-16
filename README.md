# nano-snakemake

A [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for nanopore whole genome sequencing structural variant analysis.

Contributions welcome!

## Dependencies

 - [ngmlr](https://github.com/philres/ngmlr)
 - [sniffles](https://github.com/fritzsedlazeck/Sniffles)
 - [LAST](https://github.com/mcfrith/last-rna/blob/master/last-long-reads.md) and [last-train](http://last.cbrc.jp/doc/last-train.html)
 - [minimap2](https://github.com/lh3/minimap2)
 - [NanoSV](https://github.com/mroosmalen/nanosv)
 - [svim](https://github.com/eldariont/svim)
 - [pbsv](https://github.com/PacificBiosciences/pbsv)
 - [npInv](https://github.com/haojingshao/npInv)
 - [survivor](https://github.com/fritzsedlazeck/SURVIVOR)
 - [mosdepth](https://github.com/brentp/mosdepth)
 - [cyvcf2](https://github.com/brentp/cyvcf2)
 - [vcfanno](https://github.com/brentp/vcfanno)
 - [vcftools](https://vcftools.github.io/index.html)
 - [samtools](https://github.com/samtools/samtools)
 - [matplotlib](https://github.com/matplotlib/matplotlib)
 - [seaborn](https://github.com/mwaskom/seaborn)

The required dependencies can be installed using `conda create -f environment.yaml`

## Commands
 - 'fast': minimap2 alignment with Sniffles and SVIM SV calling
 - 'precise': ngmlr alignment with Sniffles SV calling
 - 'minimap2': minimap2 alignment with Sniffles, SVIM, NanoSV and npInv SV calling
 - 'minimap2_pbsv': minimap2 alignment with pbsv-specific parameters with pbsv, SVIM, NanoSV and npInv SV calling
 - 'ngmlr': ngmlr with Sniffles, NanoSV, SVIM and npInv SV calling
 - 'last-prepare': create a LAST index and train aligner parameters using last-train


The folder "scripts" contains scripts necessary for the pipeline.   
The folder "extra_scripts" contains scripts I wrote while playing around with the results, mainly visualizations.

## Companion script
[surpyvor](https://github.com/wdecoster/surpyvor): a python wrapper around SURVIVOR, with additional convenience functions for creating high sensitivity and high confidence variant sets, calculating precision and recall metrics and visualizations using venn and upset plots.

## Citation
[Structural variants identified by Oxford Nanopore PromethION sequencing of the human genome](https://www.biorxiv.org/content/10.1101/434118v2)
