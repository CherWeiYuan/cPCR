#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 11:31:29 2025

@author: cwy
"""

from colorama import Fore, Style
from pysam import FastxFile

def parse_genome(genome_fasta):
    """
    If chrom equals a chromosome in the genome fasta, this function loads the
    fasta of chromosome in genome. If the fasta for the specific chromosome
    is not found, then this function will create all the chromosome-specific
    fasta files from the genome fasta provided. Then a genome dict will be
    loaded for the chromosome fasta.

    For analysis like minigene design, chrom can be None and the entire input
    fasta will be loaded into genome_dict.

    Input
        genome_fasta: directory and name of fasta file
        chrom       : name of chromosome
    Output
        genome_dict: dictionary of key (chromosome) and value (chromosomal
                     sequence)
    Note
        chrom.description is used instead of chrom.id to get full fasta header
        Otherwise, spaces in fasta header will lead to its truncation
    """
    genome_dict = {}

    with FastxFile(genome_fasta) as fasta_handler:
        for contig in fasta_handler:
            if contig.name not in genome_dict:
                genome_dict[str(contig.name)] = str(contig.sequence)
            else:
                raise NonUniqueFastaHeader()

    if genome_dict:
        pass
    else:
        raise SeqNotFoundError()

    return genome_dict

class NonUniqueFastaHeader(Exception):
    """
    Error for multiple of the same fasta header found in fasta
    """
    def __init__(self):
        super().__init__(Fore.RED + "Error: Non-unique fasta headers found. " +
                         Style.RESET_ALL)

class SeqNotFoundError(Exception):
    """
    Error for no sequences in fasta
    """
    def __init__(self):
        super().__init__(Fore.RED + "Error: No sequences found in fasta. " +
                         Style.RESET_ALL)