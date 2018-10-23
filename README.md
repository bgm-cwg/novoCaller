# novoCaller
## About

novoCaller is a Bayesian de novo variant calling algorithm that uses information from read-level data both in the pedigree and in unrelated samples. The method was extensively tested using large trio sequencing studies, and it consistently achieved over 98% sensitivity while giving significantly more specificity than other well known methods for the same sensitivity values. \
The method works with a vcf file alone (first layer), but it can also inspect the reads in the bam file and compute probabilities based on direction-specific (forward and reverse) reads (second layer). \
We made the software keeping VCFv4.1 in mind. The software looks at the format field of each variant call location in the VCF file and reads the AD element. Then it extracts that part from every sample.

## Compile:
```
g++ novoCaller.cpp -o novoCaller
```

## Usage:

### Run first layer (C++ code) - uses only VCF:
```
./novoCaller \
-I <path to vcf file> \
-O <path to output file for layer 1> \
-T <path to file containing sample IDs of the trios, the IDs are in the order:parent1(TAB)parent2(TAB)proband> \
-X <put 1 if you want to run on X chromosome as well, 0 otherwise> \
-P <threshold on posterior probability. Calls are made if the PP is above threshold. Use a low value like 0.005 so that a large number of calls are made for the second layer> \
-E <threshold on the ExAC allele frequency, e.g. 0.0001>
```

### Run second layer (Python code) - uses BAM files (OPTIONAL):
```
python -W ignore novoCallerBAM.py \  
-I <path to the output file from previous step (the file given in -O option)> \
-U <path to a file containing paths to the bam files from unrelated samples> \
-T <path to a file containing paths to the bam files of the trio> \
-O <path to the output file for the second layer>
```

The ignore option is given to ignore log of 0 warning.

### Example command line:

```
./novoCaller -I ./all_calls.vep.vcf -O step1_out.txt -T trio_ids.txt -X 1 -P 0.005 -E 0.008 
python -W ignore novoCallerBAM.py  -I step1_out.txt -U de_novo_unrelated_bams.txt -T de_novo_case_bams.txt -O denovo_calls.txt 
```

## Output format:
denovo_calls.txt columns are the following:
Rank, chromosome, position, reference allele, alternative allele, AF (in samples), rhos, priors, PP=posterior probability, AF_unrelated, gene_name(s)

# References
Find our paper [here](https://academic.oup.com/bioinformatics/advance-article/doi/10.1093/bioinformatics/bty749/5087716)
