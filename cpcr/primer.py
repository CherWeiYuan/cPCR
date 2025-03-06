#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 11:31:29 2025

@author: cwy
"""

import primer3

def design_primer_pair(template, ok_regions, target_regions, excluded_regions,
                       amplicon_min_size, amplicon_opt_size, amplicon_max_size,
                       primer_min_size, primer_opt_size, primer_max_size,
                       primer_min_tm, primer_opt_tm, primer_max_tm,
                       primer_gc_min, primer_gc_opt, primer_gc_max,
                       max_primer_tm_diff):
    """
    Use primer3 to design primers on the given template.
    ok_regions: a string defining allowed regions for left and right primers, e.g.
                "0,region1_length,region1_length,region2_length"
                which means left primer must lie in region 1 (starting at 0) and
                right primer in region 2 (starting at region1_length).
                
                The associated value must be a semicolon-separated list of
                <left_start>,<left_length>,<right_start>,<right_length>
                
                SEQUENCE_PRIMER_PAIR_OK_REGION_LIST=\
                    100,50,300,50 ; 900,60,, ; ,,930,100

                    Specifies that there are three choices:

                    Left primer in the 50 bp region starting at position 100 
                    AND right primer in the 50 bp region starting at position 300

                    OR

                    Left primer in the 60 bp region starting at position 900 
                    (and right primer anywhere)

                    OR

                    Right primer in the 100 bp region starting at position 930 
                    (and left primer anywhere)
    target_regions: If one or more targets is specified then a legal primer 
                    pair must flank at least one of them. A target might be a 
                    simple sequence repeat site (for example a CA repeat) or a 
                    single-base-pair polymorphism, or an exon for resequencing. 
                    The value should be a space-separated list of
                    <start>,<length>
    excluded_regions: Left and Right primers and oligos may not overlap any 
                      region specified. Space-separated list of
                      <start>,<length>

    product_size_range: tuple (min, max) of amplicon size.
    """
    # Define sequence-related settings
    seq_args = {
        'SEQUENCE_TEMPLATE': template,
        'SEQUENCE_PRIMER_PAIR_OK_REGION_LIST': ok_regions,
        'SEQUENCE_TARGET': target_regions,
        'SEQUENCE_EXCLUDED_REGION': excluded_regions 
    }

    # Define global primer design settings
    global_args = {
        'PRIMER_PRODUCT_SIZE_RANGE': [(amplicon_min_size, amplicon_max_size)],
        'PRIMER_PRODUCT_OPT_SIZE': amplicon_opt_size,

        'PRIMER_MIN_SIZE' : primer_min_size,
        'PRIMER_OPT_SIZE' : primer_opt_size,
        'PRIMER_MAX_SIZE' : primer_max_size,
        
        'PRIMER_MIN_TM'   : primer_min_tm,
        'PRIMER_OPT_TM'   : primer_opt_tm,
        'PRIMER_MAX_TM'   : primer_max_tm,
        'PRIMER_PAIR_MAX_DIFF_TM': max_primer_tm_diff,
        
        'PRIMER_MIN_GC'   : primer_gc_min,
        'PRIMER_OPT_GC_PERCENT': primer_gc_opt,
        'PRIMER_MAX_GC'   : primer_gc_max,
    }

    # Design primers
    primers = primer3.bindings.design_primers(seq_args, global_args)

    return primers

def mark_primer_duplicates(primer_df):
    """
    Mark duplicated primers
    """
    primers = primer_df.forward.to_list() +\
              primer_df.reverse.to_list()
    duplicated_primers = list(
        set([x for x in primers if primers.count(x) > 1])
        )
    primer_df["forward_duplicated"] = primer_df[
        "forward"].isin(duplicated_primers).astype(int)
    primer_df["reverse_duplicated"] = primer_df[
        "reverse"].isin(duplicated_primers).astype(int)
    return primer_df
    
