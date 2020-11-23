# Ploidetect
Tumour purity, ploidy, and copy number variation from whole-genome sequence data

## Getting started

To run Ploidetect, you will need:
- a tumor alignment bam
- a normal alignment bam

Ideally you can use the Snakemake workflow [available at this repository](https://github.com/lculibrk/Ploidetect-pipeline)

## Interpreting the output

The outputs from Ploidetect are five files:

  - models.txt
  - plots.pdf
  - cna.txt
  - cna_plots.pdf
  - meta.RDS
  
models.txt describes a number of inferred models to explain the tumor purity and ploidy of your sample. Most of the columns are gobbeldygook and need to be cleaned up, but the important columns are the tp and ploidy. The rows are ordered in the order the algorithm has scored it, so the first plot/row corresponds to the highest score.

plots.pdf is a multi-page pdf, where each page corresponds to a line in models.txt. This is intended to assist in review of the data. While great pains have been made to ensure that Ploidetect is highly accurate in model selection, it (and no tool, for that matter) is not 100% accurate and you might want to take a peek to make sure the model (red peaks) seems to agree with the distribution of observed read depths (black peaks). In tumors with a high degree of heterogeneity, the model might skip some irregular peaks. 

cna.txt contains the copy number informtation. There are a number of columns contained therein:

| column | definition |
|-|-|
| chr | chromosome |
| pos | bin START |
| end | bin END |
| segment | segment number |
| corrected_depth | corrected BIN depth |
| segment_depth | mean SEGMENT depth |
| maf | beta allele frequency at bin |
| CN | fractional copy number of segment |
| state | state (corresponds to cna_plots.pdf) |
| zygosity | HET or HOM to designate predicted zygosity |

cna_plots.pdf is a multi-page pdf containing plots of the copy number profile, with one page per chromosome.  

Currently there isn't much use for meta.RDS except for debugging.

## Troubleshooting

### My plots.pdf just shows a single blob/peak and the TP may be super low

Either the case has no CNVs or the data is too noisy at that depth level. If you're sure there are CNVs in the data, modify line 25 of the Snakefile, specifically here: 
```python3 scripts/make_windows.py - 100000``` - the 100000 is the threshold for normal depth to create variable-width bins. The size of the bins depends on the germline coverage, and Ploidetect was developed for 40x normal, 80x tumor genomes. If the ratio of tumor to normal is not 2:1, adjust the threshold accordingly. So for 40x/40x you might decide to use a larger threshold of 200000 to account for the tumor genome having more noise. 

### My copy number profile is highly fragmented/oversegmented or there are no copy number variants when there should be/undersegmented

You may have a bad estimate of tumor purity being fed to the copy number caller. Review the plots.pdf and models.txt files and find a model which explains the distribution of read depth better. Slice the row of models.txt with the correct model (row n):

```sed 'n!d' models.txt > new_models.txt```

and rerun the CNV caller with the new estimates:

```Rscript scripts/ploidetect_copynumber.R -i data/{sample_id}/segmented.RDS -m new_models.txt -p /path/to/cna_plots.pdf -o /path/to/cna.txt &> /path/to/logfile.txt"```

