# PRO-seq Expression Matching
A wrapper for Ben Martin's `expression-match.R` function tailored for outputs from the NTC PRO-seq pipeline.  This function uses quantile matching—to give similar expression distributions between control and target gene lists.

## Output
The script generates three alternative lists of expression-matched genes that are:
1.  Not differentially expressed
2.  Same number of genes
3.  Same distribution of length-normalized expression levels

*It is suggested to analyze all three to assess the quality of matching.  Downstream analysis requires subsetting matrix files using the returned gene lists.*

## Requirements
- R (version ≥4.x recommended)
- Plain text file with a single column of gene IDs of interest (e.g., upregulated genes)
- Outputs from DESeq2 (results tables and heatmap count file)
- Outputs from GGA (WORKING file)

## Usage
To run the script from the command line, 

```sh
# takes only a few seconds to run with minimal resources
# Argument order must be exactly as below
Rscript /n/data1/cores/ntc/scripts/NascentTranscriptionCore/Expression-Matching/main.R gene_ensg.txt path/to/DEseqOutput path/to/GGA.WORKING_FILE.txt <window> path/to/DEseq_metadata.txt <baseline condition name>
```
### Required Arguments (in order)
| Arg# | Example                | Description                                                                                                      |
|----------|------------------------ |------------------------------------------------------------------------------------------------------------------|
| 1        | gene_ensg.txt   | 1-column .txt with ENSG# (e.g. up- or downregulated genes)                                                                     |
| 2        | DESeqOutput/       | path to DESeqOutput folder                   |
| 3        | gga_WORKING_file.tsv    | path to GGA file with 'WORKING_FILE' in the filename (currently only works with GGA annotation)     |
| 4        | window to match expression on         | window to match expression on (e.g. if 0, matches on TSS-TES, if 250, matches on +250-TES, etc)                           |
| 5        | DEseq_metadata.txt      | Tab-delimited sample manifest for DESeq2; must include a `condition` column                                      |
| 6        | baseline condition     | Must exactly match one  conditions in DESeq metadata to use to match expression.  This is normally the WT/control to avoid matching on some weird perturbation)                               |

## Additional Details
- **Unchanged gene threshold:** |fold change| < 1.3;  padj > 0.1 (hardcoded in `main.R`)
- **Gene filtering:** Only protein-coding genes ≥1kb (after any window adjustment) are considered 
- **Expression-matching:** By default, 10 quantile groups (adjustable)
- **Replicates:** 3 seeds (1, 2, 3) for repeated sampling

The `expressionMatched` can also be used manually in an interactive R session if you want more control over the matching parameters.  See the Adelman lab's expression-matching readme for those instructions. for ease of integration with bioinformatic pipelines. It supports:

## Troubleshooting
If matching fails with an error such as:
`> Error in sample.int(length(x), size, replace, prob): cannot take a sample larger than the population when 'replace = FALSE'`
reduce the number of quantile groups (`numQuantGroups`, in the script) or decrease the control set ratio.

## Customization (make your own clone/copy of the script)
- Change thresholds, quantile bin count, or random seed count by editing `main.R`.
- Adapt column mappings if using non-standard gene annotations or experiment designs.

## Credit
Ben Martin and Geoff Nelson
