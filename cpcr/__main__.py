#!/usr/bin/env python3
"""
cDNA primer design algorithm

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
from colorama import Fore, Style
import logging
from os import makedirs
import sys
from tqdm import tqdm

from gtfparse import read_gtf
from numpy import nan
import pandas as pd
from rich_argparse import RichHelpFormatter

from cpcr.fasta import parse_genome
from cpcr.primer import design_primer_pair, mark_primer_duplicates
from cpcr.primer_blast import create_blast_database, primer_blast
from cpcr.utils import init_logging

def parse_args():
    """Parse command line arguments.
    Returns Options object with command line argument values as attributes.
    Will exit the program on a command line error.
    """
    description = "Primer design algorithm"

    # Create parent parser with common arguments
    parser = ArgumentParser(description     = description,
                            formatter_class = RichHelpFormatter)
    
    # User input parameters
    parser.add_argument("-c", "--target_exons_csv",
                        type      = str,
                        default   = None,
                        help      = \
                            "Name of tab-separated file (.tsv). \n"       +\
                            "Mandatory to have at least 5 columns: \n"    +\
                            "chrom (each row with 'chr', e.g. \n"         +\
                            "chromosome 1 will be chr1), \n"              +\
                            "start, end, strand ('+'/ '-') and \n"        +\
                            "exon_type ('noncanonical'/ 'canonical'. \n"   +\
                            "Poison exons will be added the MANE GTF \n"  +\
                            "while the canonical exons will determine \n" +\
                            "which MANE exon is the target exon.")
        
    # Fasta and GTF parameters
    parser.add_argument("-f", "--genome_fasta",
                        type      = str,
                        required  = True,
                        default   =     "genome/GRCh38.p14.genome.fasta",
                        help      = "Path and name of the genomic fasta file \n" +\
                                    "for sequence extraction from GTF.")
    parser.add_argument("-g", "--gtf_file",
                        type      = str,
                        required  = True,
                        default   = "genome/gencode.v47.annotation.gtf",
                        help      = "Path and name of the GTF file for \n " +\
                                    "for sequence extraction from genome_fasta.")
    parser.add_argument("-t", "--off_target_fasta",
                        type      = str,
                        required  = True,
                        default   = "transcriptome/gencode.v47.transcripts.fasta",
                        help      = "Path and name of the fasta file \n" +\
                                    "for off-target search (e.g. transcriptome fasta).")
    # BLAST parameters
    parser.add_argument("-b", "--blast_db_name",
                        type      = str,
                        required  = True,
                        default   = None,
                        help      = "BLAST database name.")
    parser.add_argument("-d", "--blastdb_dir",
                        type      = str,
                        default   = "./blastdb",
                        help      = "BLAST database directory.")

    parser.add_argument("--filter_low_evidence_transcripts",
                        action   = "store_true",
                        help     = "Remove low evidence transcripts from \n" +\
                                   "GENCODE GTF.")
    parser.add_argument("--constitutive_threshold",
                        type     = float,
                        default  = 0.9,
                        metavar  = "[0.0-1.0]",
                        help     = "The proportion of exons across \n"         +\
                                   "transcripts having an exon to consider \n" +\
                                   "the exon as constitutively spliced-in.\n")
    parser.add_argument("-e", "--polymerase",
                        type      = str,
                        default   = "dreamtaq",
                        choices   = ["dreamtaq", "others"],
                        help      = "Polymerase type. Used to calculate Tm. \n" +\
                                    "Dreamtaq formula is When others are used,")
    
    # Primer design parameters
    parser.add_argument("--preset_parameters",
                        action   = "store_true",
                        help     = "Use preset parameters for qPCR or splicing PCR.")
    parser.add_argument("--design_mode",
                        type      = str,
                        default   = "qPCR",
                        choices   = ["qPCR", "splicingPCR"],
                        help      = "qPCR (one primer on target exon, \n"  +\
                                    "the other on a constitutive exon) \n" +\
                                    "or splicing PCR design (primers \n"   +\
                                    "on constitutive exons flanking \n"    +\
                                    "the target exon.")
    # Preset Modes
    # ------------
    # 1. qPCR [target exon is poison exon; to check poison exon presence]
    #    preset = "qPCR_putative"
    # 2. qPCR [target exon is constitutive exon; to check mRNA expression levels 
    #    after treatment]
    #    preset = "qPCR_rnaseq"
    # 3. splicing PCR [Use MANE and alternative GENCODE transcripts]
    #    preset = "sPCR_gencode"
    # 4. splicing PCR [only 1. MANE transcript and 2. transcripts with the target exons]
    #    preset = "sPCR_rnaseq"
    parser.add_argument("--amplicon_min_size",
                        type     = int,
                        default  = 70,
                        help     = "Minimum amplicon size.")
    parser.add_argument("--amplicon_opt_size",
                        type     = int,
                        default  = 150,
                        help     = "Optimum amplicon size.")
    parser.add_argument("--amplicon_max_size",
                        type     = int,
                        default  = 250,
                        help     = "Maximum amplicon size.")
    parser.add_argument("--primer_min_size",
                        type     = int,
                        default  = 18,
                        help     = "Minimum primer size.")
    parser.add_argument("--primer_opt_size",
                        type     = int,
                        default  = 20,
                        help     = "Optimum primer size.")
    parser.add_argument("--primer_max_size",
                        type     = int,
                        default  = 24,
                        help     = "Maximum primer size.")
    parser.add_argument("--primer_min_tm",
                        type     = int,
                        default  = 58,
                        help     = "Minimum primer Tm.")
    parser.add_argument("--primer_opt_tm",
                        type     = int,
                        default  = 60,
                        help     = "Optimum primer Tm.")
    parser.add_argument("--primer_max_tm",
                        type     = int,
                        default  = 62,
                        help     = "Maximum primer Tm.")
    parser.add_argument("--max_primer_tm_diff",
                        type     = float,
                        default  = 1.0,
                        help     = "Maximum difference between primer pair Tm.")
    parser.add_argument("--primer_gc_min",
                        type     = int,
                        default  = 50,
                        help     = "Minimum primer Tm.")
    parser.add_argument("--primer_gc_opt",
                        type     = int,
                        default  = 55,
                        help     = "Optimum primer Tm.")
    parser.add_argument("--primer_gc_max",
                        type     = int,
                        default  = 60,
                        help     = "Maximum primer Tm.")

    parser.add_argument("--adjust_amplicon_size",
                        action   = "store_true",
                        help     = "Adjust amplicon size if target exon \n" +\
                                   "is small.")
    parser.add_argument("--adjust_amplicon_size_threshold",
                        type     = int,
                        default  = 80,
                        help     = "Adjust amplicon size if target exon \n" +\
                                   "is small is equal or less than this \n" +\
                                   "threshold. Used if adjust_amplicon_size is True.")
    parser.add_argument("--adjust_min_amplicon_size",
                        type     = int,
                        default  = 55,
                        help     = "Adjusted min amplicon size. Used if \n" +\
                                   "adjust_amplicon_size is True.")
    parser.add_argument("--adjust_opt_amplicon_size",
                        type     = int,
                        default  = 100,
                        help     = "Adjusted optimum amplicon size. Used if \n" +\
                                   "adjust_amplicon_size is True.")
    parser.add_argument("--adjust_max_amplicon_size",
                        type     = int,
                        default  = 200,
                        help     = "Adjusted max amplicon size. Used if \n" +\
                                   "adjust_amplicon_size is True.")
    parser.add_argument("--sanger_void",
                        type     = int,
                        default  = 30,
                        help     = "Distance (bp) from target exon that \n"   +\
                                   "does not permit primer design to \n"      +\
                                   "facilitate Sanger's sequencing of the \n" +\
                                   "target exon.")

    # Primer-blast parameters
    parser.add_argument("--total_mismatches_threshold",
                        type     = int,
                        default  = 2,
                        help     = "Total mismatches allowed in off-target binding")    
    parser.add_argument("--last_n_mismatches_threshold",
                        type     = int,
                        default  = 2,
                        help     = "Total mismatches allowed in last N bases " +\
                                   "for off-target binding.")    
    parser.add_argument("--delete_xml",
                        action   = "store_true",
                        help     = "Delete BLAST XML file after parsing.")
    parser.add_argument("--threads",
                        type     = int,
                        default  = 4,
                        help     = "Number of threads used for BLAST.")
    
    # Output parameters
    parser.add_argument("-o", "--outdir",
                        type     = str,
                        default  = "out",
                        help     = "Output directory")
    parser.add_argument("-p", "--prefix",
                        type     = str,
                        default  = "out",
                        help     = "Output file prefix")
    
    return parser.parse_args()

def main():
    """
    Main function for the primer design algorithm
    """
    args = parse_args()

    # User parameters
    target_exons_csv = args.target_exons_csv
    genome_fasta    = args.genome_fasta
    gtf_file        = args.gtf_file

    off_target_fasta = args.off_target_fasta
    blast_db_name    = args.blast_db_name
    blastdb_dir      = args.blastdb_dir

    preset_parameters = args.preset_parameters
    polymerase  = args.polymerase
    design_mode = args.design_mode
    
    filter_low_evidence_transcripts = args.filter_low_evidence_transcripts
    constitutive_threshold = args.constitutive_threshold
    
    amplicon_min_size = args.amplicon_min_size
    amplicon_opt_size = args.amplicon_opt_size
    amplicon_max_size = args.amplicon_max_size
    primer_min_size = args.primer_min_size # all PCR 18
    primer_opt_size = args.primer_opt_size # all PCR 20  
    primer_max_size = args.primer_max_size # all PCR 24
    primer_min_tm   = args.primer_min_tm # qPCR 58; sPCR 52 --> sPCR 50
    primer_opt_tm   = args.primer_opt_tm # qPCR 60; sPCR 55
    primer_max_tm   = args.primer_max_tm # qPCR 62; sPCR 60
    max_primer_tm_diff = args.max_primer_tm_diff # qPCR 1; sPCR 5
    primer_gc_min = args.primer_gc_min
    primer_gc_opt = args.primer_gc_opt
    primer_gc_max = args.primer_gc_max
    
    adjust_amplicon_size = args.adjust_amplicon_size # Implemented for splicing PCR
    adjust_amplicon_size_threshold = args.adjust_amplicon_size_threshold
    adjust_min_amplicon_size = args.adjust_min_amplicon_size
    adjust_opt_amplicon_size = args.adjust_opt_amplicon_size
    adjust_max_amplicon_size = args.adjust_max_amplicon_size
    sanger_void = args.sanger_void # Implemented for sPCR
    
    total_mismatches_threshold  = args.total_mismatches_threshold
    last_n_mismatches_threshold = args.last_n_mismatches_threshold
    
    prefix = args.prefix
    outdir = args.outdir
    delete_xml = args.delete_xml
    threads = args.threads
    
    # Make output directory
    makedirs(outdir, exist_ok = True)
    makedirs(blastdb_dir, exist_ok = True)
    
    # Init logging
    init_logging(outdir, args, silent = False)
    logging.info(Fore.GREEN + "INITIALIZE PCR DESIGN" + Style.RESET_ALL)
    version = 1.00
    logging.info(Fore.RED + f"Primer Design {version}" + Style.RESET_ALL)

    if design_mode == "qPCR" and preset_parameters:
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

    elif design_mode == "splicingPCR" and preset_parameters:
        amplicon_min_size = 100
        amplicon_opt_size = 500
        amplicon_max_size = 1500

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
        primer_gc_max = 65
    
        sanger_void = 30

    else:
        logging.info("No presets selected. Using user-defined parameters.")

    # Extract target exons from CSV
    target_exon_input_df = pd.read_csv(target_exons_csv)
    target_exon_list = []
    for t_index, _ in target_exon_input_df.iterrows():
        chrom  = target_exon_input_df.loc[t_index, "chrom"]
        if "chr" not in chrom:
            sys.exit(Fore.RED + 
                     "Error: 'chr' is not found in CSV chromosome names" +
                     Style.RESET_ALL)
        start  = int(min(target_exon_input_df.loc[t_index, "start"],
                         target_exon_input_df.loc[t_index, "end"]))
        end    = int(max(target_exon_input_df.loc[t_index, "start"],
                         target_exon_input_df.loc[t_index, "end"]))
        strand = target_exon_input_df.loc[t_index, "strand"]
        if strand not in ["+", "-"]:
            sys.exit(Fore.RED + 
                     f"Error: Strand in CSV ({strand}) is not '+' or '-'" +
                     Style.RESET_ALL)
        gene_name = target_exon_input_df.loc[t_index, "gene_name"]
        exon_type = target_exon_input_df.loc[t_index, "exon_type"]
        if exon_type not in ["noncanonical", "canonical"]:
            sys.exit(Fore.RED + f"Error: Exon type in CSV ({exon_type}) " +
                     "is not 'noncanonical' or 'canonical'" + Style.RESET_ALL)
    
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
    makedirs(blastdb_dir, exist_ok = True)

    # Create BLAST database if not done
    create_blast_database(off_target_fasta, blastdb_dir, blast_db_name)

    # Load fasta
    logging.info(Fore.GREEN + "Loading genome sequence" + Style.RESET_ALL)
    genome_dict = parse_genome(genome_fasta)

    # Load GTF
    logging.info(Fore.GREEN + "Loading GTF" + Style.RESET_ALL)
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
                    transcript_df.query(
                        f"start == {exon_start} and end == {exon_end}")
                    )
                if target_exon_present:
                    constitutive_count += 1
                else:
                    pass

            # Update MANE exons dataframe
            if constitutive_count/len(transcript_id_list) >= \
                constitutive_threshold:
                mane_exons_gtf.loc[index, "constitutive"] = 1
            else:
                mane_exons_gtf.loc[index, "constitutive"] = 0

        # Account for genes that do not meet the constitutive threshold
        if mane_exons_gtf.constitutive.sum() == 0:
            logging.warning(
                Fore.YELLOW + 
                "Warning: No MANE exon meets the constitutive threshold." +
                Style.RESET_ALL)
    
        if design_mode == "qPCR":
        # Front set
        #         ----->
        # =====Target exon=====     ==========     ==========
        #                             <-----         <-----
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
                start = min(front_exons_gtf.loc[fp_index, "start"],
                            front_exons_gtf.loc[fp_index, "end"])
                end   = max(front_exons_gtf.loc[fp_index, "start"],
                            front_exons_gtf.loc[fp_index, "end"])
                new_seq = genome_dict[chrom][start - 1: end]
                
                front_seq += new_seq
                constitutive = int(front_exons_gtf.loc[fp_index, "constitutive"])
                
                # Find ok_regions
                if constitutive:
                    # Use 0-index
                    ok_regions.append(left_ok_region + 
                                      [current_pos, len(new_seq)])
                else:
                    pass
            
                current_pos += len(new_seq)
                
            # Design front primers
            if len(ok_regions) > 0:
                primers_set1 = design_primer_pair(
                    front_seq, ok_regions, "", "",
                    amplicon_min_size, amplicon_opt_size, amplicon_max_size,
                    primer_min_size, primer_opt_size, primer_max_size,
                    primer_min_tm, primer_opt_tm, primer_max_tm,
                    primer_gc_min, primer_gc_opt, primer_gc_max, max_primer_tm_diff)
                
                logging.info(Fore.BLUE + f"{gene_name} FRONT" + Style.RESET_ALL)
                logging.info("Designed primer sets:" + Fore.GREEN + 
                             f"{primers_set1['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
                if primers_set1['PRIMER_PAIR_NUM_RETURNED'] == 0:
                    logging.info(f"Primer left explanation: {primers_set1['PRIMER_LEFT_EXPLAIN']}")
                    logging.info(f"Primer right explanation: {primers_set1['PRIMER_RIGHT_EXPLAIN']}")
                    logging.info(f"Primer pair explanation: {primers_set1['PRIMER_PAIR_EXPLAIN']}\n")
                else:
                    pass
            else:
                primers_set1 = None
                logging.warning(
                    Fore.YELLOW + 
                    f"Warning: No constitutive regions selected " +
                    "for right primer design." + Style.RESET_ALL)
            
            # Reset variables
            front_exons_gtf, target_exon_seq, front_seq = None, None, None
            current_pos, left_ok_regions, right_region, ok_regions = None, None, None, None

        # Back set
        #   ---->          ---->                                      
        # ==========     ==========     =====Target exon=====
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
                start = min(back_exons_gtf.loc[bp_index, "start"],
                            back_exons_gtf.loc[bp_index, "end"])
                end   = max(back_exons_gtf.loc[bp_index, "start"],
                            back_exons_gtf.loc[bp_index, "end"])
                new_seq = genome_dict[chrom][start - 1: end]
            
                back_seq += new_seq
                constitutive = int(back_exons_gtf.loc[bp_index, "constitutive"])
            
                # Find ok_regions
                if constitutive:
                    left_ok_regions.append([current_pos, len(new_seq)])
                else:
                    pass

                current_pos += len(new_seq)

            right_region = [current_pos, len(target_exon_seq)]
            ok_regions   = []
            for left_region in left_ok_regions:
                # If left_region is an empty list, avoid creating a new
                # ok_region
                if len(left_region) > 0:
                    ok_regions.append(left_region + right_region)
                else:
                    pass
            
            # Add target exon as the last exon
            back_seq += target_exon_seq
    
            # Design back primers
            if len(ok_regions) > 0:
                primers_set2 = design_primer_pair(
                    back_seq, ok_regions, "", "",
                    amplicon_min_size, amplicon_opt_size, amplicon_max_size,
                    primer_min_size, primer_opt_size, primer_max_size,
                    primer_min_tm, primer_opt_tm, primer_max_tm,
                    primer_gc_min, primer_gc_opt, primer_gc_max, max_primer_tm_diff)
                logging.info(Fore.BLUE + f"{gene_name} BACK" + Style.RESET_ALL)
                logging.info("Designed primer sets:" + Fore.GREEN + 
                             f"{primers_set2['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
                if primers_set2['PRIMER_PAIR_NUM_RETURNED'] == 0:
                    logging.info(f"Primer left explanation: {primers_set2['PRIMER_LEFT_EXPLAIN']}")
                    logging.info(f"Primer right explanation: {primers_set2['PRIMER_RIGHT_EXPLAIN']}")
                    logging.info(f"Primer pair explanation: {primers_set2['PRIMER_PAIR_EXPLAIN']}\n")
                else:
                    pass
            else:
                primers_set2 = None
                logging.warning(
                    Fore.YELLOW + 
                    f"Warning: No constitutive regions selected " +
                    "for left primer design." + Style.RESET_ALL)

            # Reset variables
            back_exons_gtf, target_exon_seq, back_seq = None, None, None
            current_pos, left_ok_regions, right_region, ok_regions = None, None, None, None
    
            # Get primer pair candidates
            # --------------------------
            primer_pairs1 = []
            if primers_set1 is not None:
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
                            "comments": "First exon is target"
                            }))
                        count1 += 1
                    else:
                        break
            else:
                pass

            primer_pairs2 = []
            if primers_set2 is not None:
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
                            "comments": "Last exon is the target"
                            }))
                        count2 += 1
                    else:
                        break
            else:
                pass

            primer_sets = primer_pairs1 + primer_pairs2

        # Splicing PCR
        #                   ---->                                      
        # ==========     ==========     =====Target exon=====     ==========
        #                                                           <-----
        # Note: Ignore strand because primers will be designed on both strands
        elif design_mode == "splicingPCR":
            # Process MANE exon GTF with information from target exon
            if target_exon["exon_type"] == "noncanonical":
                # Create dataframe for target exon
                target_exon_df = pd.DataFrame(
                    [[target_exon["chrom"], target_exon["start"], 
                      target_exon["end"], target_exon["strand"],
                      "noncanonical", 0]],
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
                
                # No primer should be designed on non-constitutive exons
                if not constitutive and exon_type not in [
                        "noncanonical", "canonical"]:
                    excluded_regions.append([current_pos, len(new_seq)])
                else:
                    pass
                
                # No primer should be designed on exon junctions
                excluded_regions.append([current_pos - 1, 2])
                excluded_regions.append([current_pos + len(new_seq) - 1, 2])
                
                # Avoid primer design on the target exon
                if exon_type in ["noncanonical", "canonical"]:
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
                    logging.info(Fore.YELLOW + "INFO: ADJUSTED AMPLICON SIZE" +
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
    
            logging.info("\n" + Fore.BLUE + f"{gene_name}" + Style.RESET_ALL)
            logging.info("Designed primer sets: " + Fore.GREEN + 
                  f"{primer_pairs['PRIMER_PAIR_NUM_RETURNED']}" + Style.RESET_ALL)
            if primer_pairs['PRIMER_PAIR_NUM_RETURNED'] == 0:
                logging.info(f"Primer left explanation: {primer_pairs['PRIMER_LEFT_EXPLAIN']}")
                logging.info(f"Primer right explanation: {primer_pairs['PRIMER_RIGHT_EXPLAIN']}")
                logging.info(f"Primer pair explanation: {primer_pairs['PRIMER_PAIR_EXPLAIN']}\n")
    
        else:
            sys.exit("Unexpected design mode")
    
# =============================================================================
#         # Primer-BLAST-like algorithm to check for off-targets
#         # ----------------------------------------------------
#         # Off-target search
#         # Inject a primer pair prone to off-targeting
#         primer_sets.append(("GATTCCTGGGTTTAAAAGTAAA", "CCAAAAGTAAACACTTTGTG",
#                             {"left_tm"      : None,
#                              "right_tm"     : None,
#                              "product_size" : None,
#                              "comments"     : "High off-target count control"}))
# =============================================================================

        if len(primer_sets) > 0:
            primer_df_list.append(
                primer_blast(primer_sets, gene_name, target_exon,
                            total_mismatches_threshold, last_n_mismatches_threshold,
                            amplicon_max_size, polymerase,
                            blastdb_dir, blast_db_name, outdir,
                            delete_xml, threads)
                )
        else:
            pass

        pbar.update()
    
    # Compile information into dataframe and export
    if len(primer_df_list) > 0:
        primer_df = pd.concat([x for x in primer_df_list if len(x) != 0],
                            axis = 0, ignore_index = True)
        primer_df.to_csv(f"{outdir}/{prefix}.full.tsv", sep = "\t", index = False)
    else:
        logging.error(Fore.RED + 
                      "Error: No primer is designed by primer3. Please adjust " +
                      "primer design parameters to accomodate more primers." +
                      Style.RESET_ALL)
        sys.exit("Error")

    # Filter to top primer only
    top_primer_df = primer_df.sort_values(
        by = ["gene_name", "chrom", "start", "end", 
              "off_target_count", "num_amplicons", "product_size"],
        ascending = True
        )
    top_primer_df.drop_duplicates(
        subset  = ["gene_name", "chrom", "start", "end", "strand", "comments"],
        keep    = "first", 
        inplace = True, 
        ignore_index = False
        )
    
    # Mark duplicated primers
    top_primer_df = mark_primer_duplicates(top_primer_df)
    
    # Export top primers
    top_primer_df.to_csv(f"{outdir}/{prefix}.top.tsv", sep = "\t", index = False)
    
    logging.info("\nComplete")
