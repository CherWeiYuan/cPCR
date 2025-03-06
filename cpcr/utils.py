#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar  3 11:31:29 2025

@author: cwy
"""

import logging
import sys

def init_logging(out, args, silent=False):
    """
    Initialise the logging facility. Write log statement indicating the
    program has started, and also write out the command line from sys.argv
    and allow other logging statements throughout the program to be put in
    the same file
    """
    logging.root.handlers = []

    if silent:
        logging.basicConfig(level    = logging.INFO,
                            format   = "%(asctime)s %(levelname)s - %(message)s",
                            datefmt  = "%Y-%m-%dT%H:%M:%S%z",
                            handlers = [logging.FileHandler(f"{out}/log.txt", "w+")])
    else:
        logging.basicConfig(level    = logging.INFO,
                            format   = "%(asctime)s %(levelname)s - %(message)s",
                            datefmt  = "%Y-%m-%dT%H:%M:%S%z",
                            handlers = [logging.FileHandler(f"{out}/log.txt", "w+"),
                                        logging.StreamHandler(sys.stdout)])

    logging.getLogger("tensorflow").setLevel(logging.ERROR)
    logging.info(f"Input arguments: {args}")
    #logging.info("Command line: %s", " ".join(sys.argv))
