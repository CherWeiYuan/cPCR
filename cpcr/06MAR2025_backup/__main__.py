#!/usr/bin/env python3
"""
Primer design algorithm for detecting an exon in cDNA by qPCR.
Conditions:
  - Use primer3 to design two sets of primer pairs:
      * Set 1: Forward primer binds to the target exon; Reverse primer binds to a constitutive exon.
      * Set 2: Reverse primer binds to the target exon; Forward primer binds to the constitutive exon.
  - Amplicon Size: 30–160 bp.
  - Conduct an off-target search on each primer candidate.
  
Input:
  - Target exon coordinates: dictionary with keys "chrom", "start", "end", "strand"
  - GTF file (gene models)
  - Genome FASTA file (for sequence retrieval)
  
Terminologies
  - Target exon: User-defined targeted exon, likely the poison exon
  - Primer design-allowed exon: Primers can be designed on this exon,
                                likely because it is constitutive
  - Front set: In qPCR, front set means the target exon is the first exon in 
               the template sequence
  - Back set: In qPCR, back set means the target exon is the last exon in 
               the template sequence
"""

from argparse import ArgumentParser
from argparse import RawTextHelpFormatter
from colorama import init, Fore, Style
import logging
from os import makedirs
import sys
from tqdm import tqdm

from gtfparse import read_gtf
from numpy import nan
import pandas as pd
from pcr.primer_blast import create_blast_database, primer_blast
from pcr.fasta import parse_genome
from pcr.primer import design_primer_pair

init(strip = False, convert = True, autoreset = True)
print(Fore.GREEN + "INITIALIZE PCR DESIGN" + Style.RESET_ALL)

# User parameters
target_exon_csv = "trial/input/target_exons.csv"

genome_fasta = "genome/GRCh38.p14.genome.fasta"
gtf_file = "genome/gencode.v47.annotation.gtf"

off_target_fasta = "transcriptome/gencode.v47.transcripts.fasta" # "transcriptome/gencode.v47.pc_transcripts.fasta"
blast_db_name = "gencode.v47.transcripts" # "gencode.v47.pc_transcripts"
blastdb_dir = "./blastdb"

prefix = "trial"
outdir = "./trial"

preset = "sPCR_gencode" # qPCR_putative, sPCR_gencode, qPCR_rna_seq
polymerase = "dreamtaq"
design_mode = "qPCR"

filter_low_evidence_transcripts = False
constitutive_threshold = 0.8

amplicon_min_size = 70
amplicon_opt_size = 150
amplicon_max_size = 250

primer_min_size = 18 # all PCR 18
primer_opt_size = 20 # all PCR 20       
primer_max_size = 24 # all PCR 24

primer_min_tm   = 58 # qPCR 58; sPCR 52 --> sPCR 50
primer_opt_tm   = 60 # qPCR 60; sPCR 55
primer_max_tm   = 62 # qPCR 62; sPCR 60

max_primer_tm_diff = 1 # qPCR 1; sPCR 5

primer_gc_min = 50
primer_gc_opt = 55
primer_gc_max = 60

adjust_amplicon_size = False # Implemented for splicing PCR
adjust_amplicon_size_threshold = 80
adjust_min_amplicon_size = 100
adjust_opt_amplicon_size = 100
adjust_max_amplicon_size = 200

sanger_void = 0 # Implemented for sPCR

total_mismatches_threshold  = 2
last_n_mismatches_threshold = 5

delete_xml = True
threads = 8

#### 4 Modes
# 1. qPCR [target exon is poison exon; to check poison exon presence]
#    preset = "qPCR_putative"
# 2. qPCR [target exon is constitutive exon; to check mRNA expression levels 
#    after treatment]
#    preset = "qPCR_rnaseq"
# 3. splicing PCR [Use MANE and alternative GENCODE transcripts]
#    preset = "sPCR_gencode"
# 4. splicing PCR [only 1. MANE transcript and 2. transcripts with the target exons]
#    preset = "sPCR_rnaseq"
if preset == "qPCR_putative" or preset == "qPCR_rna_seq":
    design_mode = "qPCR"
    amplicon_min_size = 70
    amplicon_opt_size = 150
    amplicon_max_size = 250

    adjust_amplicon_size = False

    primer_min_size = 18
    primer_opt_size = 20      
    primer_max_size = 24

    primer_min_tm = 58
    primer_opt_tm = 60
    primer_max_tm = 62
    max_primer_tm_diff = 1

    primer_gc_min = 50
    primer_gc_opt = 55
    primer_gc_max = 60
    sanger_void = 0

elif preset == "sPCR_gencode":
    design_mode = "splicing_pcr"
    amplicon_min_size = 100
    amplicon_opt_size = 500
    amplicon_max_size = 2000

    adjust_amplicon_size = True
    adjust_amplicon_size_threshold = 80
    adjust_min_amplicon_size = 100
    adjust_opt_amplicon_size = 200
    adjust_max_amplicon_size = 300

    primer_min_size = 18
    primer_opt_size = 20      
    primer_max_size = 24

    primer_min_tm = 50
    primer_opt_tm = 55
    primer_max_tm = 60
    max_primer_tm_diff = 5

    primer_gc_min = 40
    primer_gc_opt = 50
    primer_gc_max = 60

    sanger_void = 30

else:
    print("No presets selected. Using user-defined parameters.")

# Extract target exons from CSV
target_exon_input_df = pd.read_csv(target_exon_csv)
target_exon_list = []
for t_index, _ in target_exon_input_df.iterrows():
    chrom  = target_exon_input_df.loc[t_index, "chrom"]
    if "chr" not in chrom:
        sys.exit(Fore.RED + 
                 "Error: 'chr' is not found in CSV chromosome names" +
                 Style.RESET_ALL)
    start  = int(target_exon_input_df.loc[t_index, "start"])
    end    = int(target_exon_input_df.loc[t_index, "end"])
    strand = target_exon_input_df.loc[t_index, "strand"]
    if strand not in ["+", "-"]:
        sys.exit(Fore.RED + 
                 f"Error: Strand in CSV ({strand}) is not '+' or '-'" +
                 Style.RESET_ALL)
    gene_name = target_exon_input_df.loc[t_index, "gene_name"]
    exon_type = target_exon_input_df.loc[t_index, "exon_type"]
    if exon_type not in ["poison_exon", "canonical"]:
        sys.exit(Fore.RED + f"Error: Exon type in CSV ({exon_type}) " +
                 "is not 'poison_exon' or 'canonical'" + Style.RESET_ALL)

    target_exon_list.append(
        {"gene_name" : gene_name, 
         "chrom"     : chrom, 
         "start"     : start, 
         "end"       : end, 
         "strand"    : strand,
         "exon_type" : exon_type}
        )

# Make output directory
makedirs(outdir, exist_ok = True)
makedirs(f"{outdir}/pcr_temp", exist_ok = True)
makedirs(blastdb_dir, exist_ok = True)

# Create BLAST database if not done
create_blast_database(off_target_fasta, blastdb_dir, blast_db_name)

# Load fasta
print(Fore.GREEN + "Loading genome sequence" + Style.RESET_ALL)
genome_dict = parse_genome(genome_fasta)

# Load GTF
print(Fore.GREEN + "Loading GTF" + Style.RESET_ALL)
gtf = read_gtf(gtf_file).to_pandas()

# Filter GTF to keep only transcripts and exons
# Avoid design in UTR because of UTR length uncertainty
gtf = gtf[gtf.feature.isin(["exon"])]

if filter_low_evidence_transcripts:
    # Remove transcripts with weak support level
    gtf = gtf.query("level == '1' or level == '2'")
    gtf = gtf.query("transcript_support_level == '1' or transcript_support_level == '2'")
else:
    pass

# Start loop per target exon
primer_df_list = []
pbar = tqdm(total = len(target_exon_list), position = 0, leave = True)
for target_exon in target_exon_list:

    # Subset GTF to gene name
    gene_name = target_exon["gene_name"]
    gene_gtf = gtf.query(f"gene_name == '{gene_name}'")
    
    # Subset GTF to MANE
    # Design primers only on the MANE transcript
    # Why? If an exon is not on MANE, it is unlikely to be constitutive anyway
    mane_gtf = gene_gtf[gene_gtf.tag.str.contains("MANE_Select")]
    mane_exons_gtf = mane_gtf.query("feature == 'exon'")
    
    # Sort by increasing genomic coordinates
    mane_exons_gtf.sort_values(
        by = ["seqname", "start", "end"], 
        ascending = True, 
        inplace   = True,
        ignore_index = True)

    # Get a list of gene transcript ID
    transcript_id_list = gene_gtf.transcript_id.unique().tolist()
    
    # Mark exons that are [1] MANE and [2] constitutive
    # Definition of "constitutive" is the percentage of transcripts with the
    # exon is above the set constitutive_threshold
    constitutive_count = 0
    for index, _ in mane_exons_gtf.iterrows():
        exon_start = int(mane_exons_gtf.loc[index, "start"])
        exon_end   = int(mane_exons_gtf.loc[index, "end"])
    
        # Assume exon is constitutive and attempt to disprove it
        for tid in transcript_id_list:
            transcript_df = gene_gtf.query(f"transcript_id == '{tid}'")
            target_exon_present = len(
                transcript_df.query(f"start == {exon_start} and end == {exon_end}")
                )
            if target_exon_present:
                constitutive_count += 1
            else:
                pass
        
        # Update MANE exons dataframe
        if constitutive_count/len(transcript_id_list) >= constitutive_threshold:
            mane_exons_gtf.loc[index, "constitutive"] = 1
        else:
            mane_exons_gtf.loc[index, "constitutive"] = 0
    
    # Account for genes that do not meet the constitutive threshold
    if mane_exons_gtf.constitutive.sum() == 0:
        sys.exit(Fore.RED + 
                 "Error: No MANE exon meets the constitutive threshold" +
                 Style.RESET_ALL)

    if design_mode == "qPCR":
        # Front set
        #                        ----->
        # ==========     =====Target exon=====     ==========     ==========
        #                                            <-----         <-----
        # Note: Ignore strand because primers will be designed on both strands
        
        # Assemble sequence where the target exon is the first exon
        front_exons_gtf = mane_exons_gtf.query(f"start > {target_exon['end']}")
        
        front_seq = genome_dict[target_exon['chrom']][
            target_exon['start'] - 1: target_exon['end']]
        
        current_pos = len(front_seq) - 1 # Use 0-index
        left_ok_region = [0, len(front_seq)]
        ok_regions = []
        for fp_index, _ in front_exons_gtf.iterrows():
            chrom = front_exons_gtf.loc[fp_index, "seqname"]
            start = front_exons_gtf.loc[fp_index, "start"]
            end   = front_exons_gtf.loc[fp_index, "end"]
            new_seq = genome_dict[chrom][start - 1: end]
            
            front_seq += new_seq
            constitutive = front_exons_gtf.loc[fp_index, "constitutive"]
            
            # Find ok_regions
            if constitutive:
                # Use 0-index
                ok_regions.append(left_ok_region + 
                                  [current_pos, len(new_seq)])
            else:
                pass
        
            current_pos += len(new_seq)
            
        # Design front primers
        primers_set1 = design_primer_pair(
            front_seq, ok_regions, "", "",
            amplicon_min_size, amplicon_opt_size, amplicon_max_size,
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            primer_gc_min, primer_gc_opt, primer_gc_max, max_primer_tm_diff)
        
        # Reset variables
        front_exons_gtf, target_exon_seq, front_seq = None, None, None
        current_pos, left_ok_regions, right_region, ok_regions = None, None, None, None
        
        print(Fore.BLUE + f"{gene_name} FRONT" + Style.RESET_ALL)
        print("\n" + Fore.BLUE + len(f"{gene_name} FRONT") * "-" + Style.RESET_ALL)
        print("Designed primer sets:" + Fore.GREEN + 
              f"{primers_set1['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
        if primers_set1['PRIMER_PAIR_NUM_RETURNED'] == 0:
            print(f"Primer left explanation: {primers_set1['PRIMER_LEFT_EXPLAIN']}")
            print(f"Primer right explanation: {primers_set1['PRIMER_RIGHT_EXPLAIN']}")
            print(f"Primer pair explanation: {primers_set1['PRIMER_PAIR_EXPLAIN']}\n")
        else:
            pass

        # Back set
        #   ---->          ---->                                      
        # ==========     ==========     =====Target exon=====     ==========
        #                                      <-----
        # Note: Ignore strand because primers will be designed on both strands
        
        # Assemble sequence where the target exon is the last exon
        back_exons_gtf = mane_exons_gtf.query(f"end < {target_exon['start']}")
        
        # Get target exon sequence
        target_exon_seq = genome_dict[target_exon['chrom']][
            target_exon['start'] - 1: target_exon['end']]
        
        back_seq = ""
        
        current_pos = 0
        left_ok_regions = []
        for bp_index, _ in back_exons_gtf.iterrows():
            chrom = back_exons_gtf.loc[bp_index, "seqname"]
            start = back_exons_gtf.loc[bp_index, "start"]
            end   = back_exons_gtf.loc[bp_index, "end"]
            new_seq = genome_dict[chrom][start - 1: end]
        
            back_seq += new_seq
            constitutive = back_exons_gtf.loc[bp_index, "constitutive"]
        
            # Find ok_regions
            if constitutive:
                left_ok_regions.append([current_pos, len(new_seq)])
            else:
                pass
        
            current_pos += len(new_seq)
    
        right_region = [current_pos, len(target_exon_seq)]
        ok_regions   = []
        for left_region in left_ok_regions:
            ok_regions.append(left_region + right_region)
    
        # Add target exon as the last exon
        back_seq += target_exon_seq

        # Design back primers
        primers_set2 = design_primer_pair(
            back_seq, ok_regions, "", "",
            amplicon_min_size, amplicon_opt_size, amplicon_max_size,
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            primer_gc_min, primer_gc_opt, primer_gc_max, max_primer_tm_diff)
        
        # Reset variables
        back_exons_gtf, target_exon_seq, back_seq = None, None, None
        current_pos, left_ok_regions, right_region, ok_regions = None, None, None, None

        print(Fore.BLUE + f"{gene_name} BACK" + Style.RESET_ALL)
        print("\n" + Fore.BLUE + len(f"{gene_name} BACK") * "-" + Style.RESET_ALL)
        print("Designed primer sets:" + Fore.GREEN + 
              f"{primers_set2['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
        if primers_set2['PRIMER_PAIR_NUM_RETURNED'] == 0:
            print(f"Primer left explanation: {primers_set2['PRIMER_LEFT_EXPLAIN']}")
            print(f"Primer right explanation: {primers_set2['PRIMER_RIGHT_EXPLAIN']}")
            print(f"Primer pair explanation: {primers_set2['PRIMER_PAIR_EXPLAIN']}\n")
        else:
            pass

        # Get primer pair candidates
        # --------------------------
        primer_pairs1 = []
        count1 = 0
        while True:
            if f"PRIMER_LEFT_{count1}_SEQUENCE" in primers_set1 and \
               f"PRIMER_RIGHT_{count1}_SEQUENCE" in primers_set1:
                primer_pairs1.append(
                    (primers_set1[f"PRIMER_LEFT_{count1}_SEQUENCE"],
                     primers_set1[f"PRIMER_RIGHT_{count1}_SEQUENCE"],
                     {"left_tm": primers_set1[f"PRIMER_LEFT_{count1}_TM"],
                      "right_tm": primers_set1[f"PRIMER_RIGHT_{count1}_TM"],
                      "product_size": primers_set1[f"PRIMER_PAIR_{count1}_PRODUCT_SIZE"],
                      "comments": "Front exon set"
                      }))
                count1 += 1
            else:
                break
            
        primer_pairs2 = []
        count2 = 0
        while True:
            if f"PRIMER_LEFT_{count2}_SEQUENCE" in primers_set2 and \
               f"PRIMER_RIGHT_{count2}_SEQUENCE" in primers_set2:
                primer_pairs2.append(
                    (primers_set2[f"PRIMER_LEFT_{count2}_SEQUENCE"],
                     primers_set2[f"PRIMER_RIGHT_{count2}_SEQUENCE"],
                     {"left_tm": primers_set2[f"PRIMER_LEFT_{count2}_TM"],
                      "right_tm": primers_set2[f"PRIMER_RIGHT_{count2}_TM"],
                      "product_size": primers_set2[f"PRIMER_PAIR_{count2}_PRODUCT_SIZE"],
                      "comments": "Back exon set"
                      }))
                count2 += 1
            else:
                break

        primer_sets = primer_pairs1 + primer_pairs2

    elif design_mode == "splicing_pcr":
        # Process MANE exon GTF with information from target exon
        if target_exon["exon_type"] == "poison_exon":
            # Create dataframe for target exon
            target_exon_df = pd.DataFrame(
                [[target_exon["chrom"], target_exon["start"], 
                  target_exon["end"], target_exon["strand"],
                  "poison_exon", 0]],
                columns = ["seqname", "start", "end", "strand", "EXON_TYPE",
                           "constitutive"]
                )

            # Concatenate with MANE GTF
            mane_exons_gtf = pd.concat([mane_exons_gtf, target_exon_df], axis = 0)
            mane_exons_gtf.sort_values(by = ["seqname", "start", "end"],
                                       inplace = True, ignore_index = True)

        elif target_exon["exon_type"] == "canonical":
            # Check if target exon intersects mane exons
            # If yes, make the mane exons canonical
            start1, end1 = target_exon["start"], target_exon["end"]
            for m_index, _ in mane_exons_gtf.iterrows():
                start2 = mane_exons_gtf.loc[m_index, "start"]
                end2   = mane_exons_gtf.loc[m_index, "end"]
                if max(start1, start2) < min(end1, end2):
                    mane_exons_gtf.loc[m_index, "EXON_TYPE"] = "canonical"

        else:
            sys.exit(Fore.RED + "Error: unexpected exon type" +
                     Style.RESET_ALL)
            
        # Extract included/ excluded regions for primer3 design
        current_pos = 0
        mane_seq    = ""
        min_target_exon_size = None
        excluded_regions = []
        primer3_sequence_target = []
        for m_index, _ in mane_exons_gtf.iterrows():
            chrom = mane_exons_gtf.loc[m_index, "seqname"]
            start = mane_exons_gtf.loc[m_index, "start"]
            end   = mane_exons_gtf.loc[m_index, "end"]
            exon_type = mane_exons_gtf.loc[m_index, "EXON_TYPE"]
            constitutive = mane_exons_gtf.loc[m_index, "constitutive"]
            
            new_seq = genome_dict[chrom][start - 1: end]
            mane_seq += new_seq
            
            if not constitutive and exon_type not in [
                    "poison_exon", "canonical"]:
                excluded_regions.append([current_pos, len(new_seq)])
            else:
                pass

            if exon_type in ["poison_exon", "canonical"]:
                primer3_sequence_target.append([current_pos, len(new_seq)])
                
                # Leave space before and after target exons for Sanger's seq
                excluded_regions.append(
                    [current_pos - sanger_void, current_pos])
                
                excluded_regions.append(
                    [current_pos + len(new_seq), 
                     current_pos + len(new_seq) + sanger_void])

                # Get the size of the smallest target exon
                if min_target_exon_size is None:
                    min_target_exon_size = len(new_seq)
                elif min_target_exon_size > len(new_seq):
                    min_target_exon_size = len(new_seq)
                else:
                    pass
            else:
                pass
            
            current_pos += len(new_seq)

        # Remove excluded regions that extend beyond sequence length
        cleaned_excluded_regions = []
        for region in excluded_regions:
            current_pos = region[0]
            length      = region[1]
            if current_pos + length <= len(mane_seq):
                cleaned_excluded_regions.append(region)
            else:
                pass
        excluded_regions = cleaned_excluded_regions
        
        # Check if target sequence is found
        if not primer3_sequence_target:
            sys.exit(Fore.RED +"Error: Target exon not found" +
                     Style.RESET_ALL)
            
        # Decrease amplicon_max_size if target exon is small
        # This is to allow gel bands with small size difference
        # to be differentiated
        if adjust_amplicon_size and \
            min_target_exon_size <= adjust_amplicon_size_threshold:
                spcr_amp_min_size = adjust_min_amplicon_size
                spcr_amp_opt_size = adjust_opt_amplicon_size
                spcr_amp_max_size = adjust_max_amplicon_size
                print(Fore.YELLOW + "INFO: ADJUSTED AMPLICON SIZE" +
                      Style.RESET_ALL)
        else:
            spcr_amp_min_size = amplicon_min_size
            spcr_amp_opt_size = amplicon_opt_size
            spcr_amp_max_size = amplicon_max_size

        # Design primers
        primer_pairs = design_primer_pair(
            mane_seq, "", 
            primer3_sequence_target, 
            excluded_regions,
            spcr_amp_min_size, spcr_amp_opt_size, spcr_amp_max_size,
            primer_min_size, primer_opt_size, primer_max_size,
            primer_min_tm, primer_opt_tm, primer_max_tm,
            primer_gc_min, primer_gc_opt, primer_gc_max, max_primer_tm_diff)

        # Get primer sets
        primer_sets = []
        count = 0
        while True:
            if f"PRIMER_LEFT_{count}_SEQUENCE" in primer_pairs and \
               f"PRIMER_RIGHT_{count}_SEQUENCE" in primer_pairs:
                primer_sets.append(
                    (primer_pairs[f"PRIMER_LEFT_{count}_SEQUENCE"],
                     primer_pairs[f"PRIMER_RIGHT_{count}_SEQUENCE"],
                     {"left_tm": primer_pairs[f"PRIMER_LEFT_{count}_TM"],
                      "right_tm": primer_pairs[f"PRIMER_RIGHT_{count}_TM"],
                      "product_size": primer_pairs[f"PRIMER_PAIR_{count}_PRODUCT_SIZE"],
                      "comments": ""
                      }))
                count += 1
            else:
                break

        print("\n" + Fore.BLUE + f"{gene_name}" + Style.RESET_ALL)
        print("\n" + Fore.BLUE + len(gene_name) * "-" + Style.RESET_ALL)
        print("Designed primer sets: " + Fore.GREEN + 
              f"{primer_pairs['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
        if primer_pairs['PRIMER_PAIR_NUM_RETURNED'] == 0:
            print(f"Primer left explanation: {primer_pairs['PRIMER_LEFT_EXPLAIN']}")
            print(f"Primer right explanation: {primer_pairs['PRIMER_RIGHT_EXPLAIN']}")
            print(f"Primer pair explanation: {primer_pairs['PRIMER_PAIR_EXPLAIN']}\n")

    else:
        sys.exit("Unexpected design mode")

    # Primer-BLAST-like algorithm to check for off-targets
    # ----------------------------------------------------
    # Off-target search
    # Inject a primer pair prone to off-targeting
    primer_sets.append(("GATTCCTGGGTTTAAAAGTAAA", "CCAAAAGTAAACACTTTGTG",
                        {"left_tm"      : None,
                         "right_tm"     : None,
                         "product_size" : None,
                         "comments"     : "High off-target count control"}))

    primer_df_list.append(
        primer_blast(primer_sets, gene_name, target_exon,
                     total_mismatches_threshold, last_n_mismatches_threshold,
                     amplicon_max_size, polymerase,
                     blastdb_dir, blast_db_name, outdir,
                     delete_xml, threads)
        )
    
    pbar.update()

# Compile information into dataframe and export
primer_df = pd.concat(primer_df_list, axis = 0, ignore_index = True)
primer_df.to_csv(f"{outdir}/{prefix}.full.tsv", sep = "\t", index = False)

# Filter to top primer only
top_primer_df = primer_df.sort_values(
    by = ["gene_name", "chrom", "start", "end", "off_target_count", "num_amplicons"],
    ascending = True
    )
top_primer_df.drop_duplicates(
    subset  = ["gene_name", "chrom", "start", "end", "strand", "comments"],
    keep    = "first", 
    inplace = True, 
    ignore_index = False
    )
top_primer_df.to_csv(f"{outdir}/{prefix}.top.tsv", sep = "\t", index = False)

print("\nComplete")
