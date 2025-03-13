# cPCR
cDNA PCR primer design

# Set up
## Create mamba environment
```
mamba create --name cpcr \
    blast biopython gtfparse rich-argparse spyder primer3 primer3-py pysam tqdm \
    -c bioconda -y
```

## Install cPCR
```
git clone https://github.com/CherWeiYuan/Panthera.Git
pip install -e .
```

# How it works?
The algorithm designs PCR primers on constitutive exons (defined as at least 80% of the GTF's transcript having the exon or otherwise specified by --constitutive_threshold) using the genome fasta and genome GTF specified by the user (recommended to use GENCODE files) via the primer3 algorithm.

A set of BLAST functions is written to check the primers designed for off-target binding and the creation of rogue amplicons. Primer off-target binding is defined as the default in NCBI's Primer-BLAST: "Primer must have at least 2 total mismatches to unintended targets, including at least 2 mismatches within the last 5 bps at the 3' end.  Ignore targets that have 6 or more mismatches to the primer."

The algorithm works in two modes: qPCR and splicing PCR. The former, qPCR, design primers on the target exon and a flanking constitutive exon. The latter, splicing PCR, design primers on exons flanking the target exon.

# How to run?
The parameters of the algorithm can be found with cpcr --help.

You must specify the following:

[1] --target_exons_csv: A comma-separated file with 5 columns: chrom (name based on the GTF), start, end, strand (-/+) and exon_type (canonical/ noncanonical).
Each row will represent a target exon's information. See "/test_data/csv/target_exons.csv" for an example.
If the exon_type is canonical, the algorithm finds the exon(s) in the MANE transcript and mark it as a target exon.
If the exon_type is noncanonical, the algorithm adds the target exon to the MANE transcript.

[2] --design_mode: If the method is qPCR, one primer will be designed on the target exon and the other on a constitutive exon.
If the method is splicing_pcr, then the primers will be designed on exons that flank the target exon. The exon with the primers may not be immediately adjacent to the target exon.

Default parameters run:
```
mamba activate cpcr
cpcr \
  --target_exons_csv [csv file] \
  --design_mode [qPCR/ splicingPCR] \
  --preset_parameters \
  --outdir [output directory] \
  --genome_fasta ./genome/GRCh38.p14.genome.fasta \
  --gtf_file ./genome/gencode.v47.annotation.gtf \
  --off_target_fasta ./transcriptome/gencode.v47.transcripts.fasta \
  --blast_db_name gencode.v47.transcripts \
  --delete_xml \
  --threads 4
mamba deactivate
```

# Output
The output will be deposited in the folder specified in the --outdir. Two TSVs will be created.
The TSV with suffix ".full.tsv" contains all the primers designed while the suffix ".top.tsv" contain the top primers.
For examples, see "/test_data/output_qPCR" for qPCR and "/test_data/output_splicingPCR" for splicing PCR.
