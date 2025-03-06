#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  4 14:41:54 2025

@author: cwy
"""

from os import path, remove

from Bio.Blast.Applications import NcbimakeblastdbCommandline, NcbiblastnCommandline
from Bio.Blast import NCBIXML
from numpy import nan
import pandas as pd

def create_blast_database(genome_fasta, blastdb_dir, blast_db_name):
    """
    Create BLAST database
    """
    if not path.isfile(f"{blastdb_dir}/{blast_db_name}.nsq"):
        print("Creating BLAST database.")
        makeblastdb_cline = NcbimakeblastdbCommandline(
            input_file = genome_fasta,
            dbtype     = "nucl",
            out        = f"{blastdb_dir}/{blast_db_name}")
        makeblastdb_cline()

def calculate_mismatches(hsp):
    """
    Count number of mismatches in high-scoring pairs
    """
    mismatches = 0
    for q, s in zip(hsp.query, hsp.sbjct):
        if q != s or q == '-' or s == '-':
            mismatches += 1
        else:
            pass
    return mismatches

def calculate_last_n_mismatches(hsp, n = 5):
    """
    Count number of mismatches in the last N nucleotides of the alignment
    """
    query_seq    = hsp.query.replace("-", "")
    last_n_query = query_seq[-n:] if len(query_seq) >= n else query_seq
    last_n_sbjct = hsp.sbjct[-len(last_n_query):]

    mismatches = 0
    for q, s in zip(last_n_query, last_n_sbjct):
        if q != s or q == '-' or s == '-':
            mismatches += 1

    return mismatches

def run_blast(primer_seq, primer_name, strand, blast_db_name, threads, outdir):
    """
    Run local BLAST
    """
    xml = f"{outdir}/pcr_temp/{primer_name}.blast.xml"
    cline = NcbiblastnCommandline(
        db     = blast_db_name,
        out    = xml,
        outfmt = 5,
        word_size = 7,
        task      = "blastn-short",
        evalue    = 10,
        strand    = strand,
        dust         = "'no'",    # No masking of any seqs because any off-
        soft_masking = "'false'", # target will result in the lost of depth
        num_alignments = 1000,
        num_threads    = threads)

    # Run BLAST which output XML
    stdout, stderr = cline(stdin = primer_seq)
    
    return xml

def parse_blast(xml,
                total_mismatches_threshold  = 2,
                last_n_mismatches_threshold = 5,
                delete_xml = True):
    """
    Parse BLAST output in XML file
    """
    hits = []
    result_handle = open(xml)
    blast_records = NCBIXML.read(result_handle)
    for alignment in blast_records.alignments:
        title = alignment.title
        for hsp in alignment.hsps:
            # Calculate mismatches
            total_mismatches = calculate_mismatches(hsp)
            if total_mismatches >= 6:
                continue
            last_n_mismatches = calculate_last_n_mismatches(hsp)
            if total_mismatches <= total_mismatches_threshold and \
                last_n_mismatches <= last_n_mismatches_threshold:
                    chrom = alignment.title.split()[0]
                    start = min(hsp.sbjct_start, hsp.sbjct_end)
                    end   = max(hsp.sbjct_start, hsp.sbjct_end)

                    hits.append({
                        "title" : title,
                        "chrom" : chrom,
                        "start" : start,
                        "end"   : end,
                        "total_mismatches"  : total_mismatches,
                        "last_n_mismatches" : last_n_mismatches})
    
    # Delete XML file
    if delete_xml: 
        remove(xml)

    return hits

def find_amplicons(forward_hits, reverse_hits, amplicon_max_size):
    """
    Find amplicons given forward primer, reverse primer binding sites and
    max accepted amplicon size
    """
    amplicons = []
    for f_hit in forward_hits:
        for r_hit in reverse_hits:
            if f_hit["chrom"] != r_hit["chrom"]:
                continue
            
            if f_hit["start"] < r_hit["end"] and \
               (r_hit["end"] - f_hit["start"] + 1) <= amplicon_max_size:
                
                amplicons.append({
                    "title"  : f_hit["title"],
                    "chrom"  : f_hit["chrom"],
                    "start"  : f_hit["start"],
                    "end"    : r_hit["end"],
                    "length" : r_hit["end"] - f_hit["start"] + 1,
                    "f_mismatches" : f_hit["total_mismatches"],
                    "r_mismatches" : r_hit["total_mismatches"]})

    return amplicons

def tm_calculator(seq, polymerase):
    """
    Calculate Tm based on polymerase type
    """
    seq = seq.upper()
    if polymerase == "dreamtaq":
        gc = 4 * (seq.count("G") + seq.count("C"))
        at = 2 * (seq.count("A") + seq.count("T"))
        return gc + at
    else:
        return nan

def primer_blast(primer_sets, gene_name, target_exon,
                 total_mismatches_threshold, last_n_mismatches_threshold,
                 amplicon_max_size, polymerase,
                 blastdb_dir, blast_db_name, outdir,
                 delete_xml, threads):
    """
    Run PRIMER-BLAST-like function to 
    [1] identify primer off-target binding and
    [2] use the binding sites to identify off-target amplicons
    """
    primer_df_seed = []
    for forward_primer, reverse_primer, info_dict in primer_sets:
        # Forward primer blast on plus strand
        forward_plus_xml  = run_blast(
            forward_primer, "forward_plus", "plus",
            f"{blastdb_dir}/{blast_db_name}", threads, outdir)
        forward_plus_hits = parse_blast(
            forward_plus_xml, 
            total_mismatches_threshold  = total_mismatches_threshold,
            last_n_mismatches_threshold = last_n_mismatches_threshold,
            delete_xml = delete_xml
            )
    
        # Reverse primer blast on minus strand
        reverse_minus_xml  = run_blast(
            reverse_primer, "reverse_minus", "minus", 
            f"{blastdb_dir}/{blast_db_name}", threads, outdir)
        reverse_minus_hits = parse_blast(
            reverse_minus_xml,
            total_mismatches_threshold  = total_mismatches_threshold,
            last_n_mismatches_threshold = last_n_mismatches_threshold,
            delete_xml = delete_xml)
    
        # Forward primer blast on minus strand
        forward_minus_xml  = run_blast(
            forward_primer, "forward_minus", "minus",
            f"{blastdb_dir}/{blast_db_name}", threads, outdir)
        forward_minus_hits = parse_blast(
            forward_minus_xml, 
            total_mismatches_threshold  = total_mismatches_threshold,
            last_n_mismatches_threshold = last_n_mismatches_threshold,
            delete_xml = delete_xml
            )
        
        # Reverse primer blast on minus strand
        reverse_plus_xml  = run_blast(
            reverse_primer, "reverse_plus", "plus", 
            f"{blastdb_dir}/{blast_db_name}", threads, outdir)
        reverse_plus_hits = parse_blast(
            reverse_plus_xml,
            total_mismatches_threshold  = total_mismatches_threshold,
            last_n_mismatches_threshold = last_n_mismatches_threshold,
            delete_xml = delete_xml)
    
        # Compute the possible amplicons from the blast hits
        amplicons = []
        amplicons += find_amplicons(forward_plus_hits, reverse_minus_hits, 
                                    amplicon_max_size)
        amplicons += find_amplicons(reverse_plus_hits, forward_minus_hits,
                                    amplicon_max_size)
        
        # Find the number of amplicons not in the target gene
        off_target_count = 0
        for product in amplicons:
            if gene_name not in product["title"]:
                off_target_count += 1
            else:
                pass

        # Save primer pair information
        primer_df_seed.append(
            (target_exon["gene_name"], target_exon["chrom"], 
             target_exon["start"], target_exon["end"], target_exon["strand"], 
             forward_primer, reverse_primer,
             tm_calculator(forward_primer, polymerase), 
             tm_calculator(reverse_primer, polymerase),
             info_dict["left_tm"], info_dict["right_tm"],
             info_dict["product_size"], off_target_count, 
             len(amplicons), amplicons, info_dict["comments"]))
    
    # Output dataframe
    primer_df = pd.DataFrame(
        primer_df_seed,
        columns = ["gene_name", "chrom", "start", "end", "strand", 
                   "forward", "reverse",
                   f"forward_{polymerase}_tm", f"reverse_{polymerase}_tm", 
                   "forward_primer3_tm", "reverse_primer3_tm", 
                   "product_size", "off_target_count", 
                   "num_amplicons", "amplicons", "comments"])
    
    return primer_df