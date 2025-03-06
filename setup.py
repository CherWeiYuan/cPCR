from setuptools import setup, find_packages

LONG_DESCRIPTION = \
"""
cDNA PCR primer design
"""

setup(name         = "cpcr",
      author       = "CHER_WEI_YUAN",
      version      = "1.0.0",
      author_email = "E0031403@U.NUS.EDU",
      packages     = find_packages(include=["cpcr"]),
      entry_points = {"console_scripts": ["cpcr = cpcr.__main__: main"]},
      description  = ("Bioinformatics commandline tool"),
      long_description = (LONG_DESCRIPTION))
