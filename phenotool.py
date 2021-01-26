#!/home/fls530/anaconda3/bin/python

# An Idea

Version = """0.4 (Development Version)
"""

# A handy tool to fix and re-afirm sample files for GWAS
#	should be able to output .sam and .psam formats (the latter being so much easier...)

# We should define a phenotype class

EPILOG_Legal = """
Copyright 2019 Carsten Friis Rundsten <fls530@ku.dk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
MA 02110-1301, USA.
"""

OPTION_files="""
Input File(s). %(prog)s accepts one or more input files including '-' symbolizing stdin. The precise format of each input file
will be autodetected, but it should be some form of delimited text data file like 'csv' or tab-delimited.
"""

OPTION_samples="""
File with samples to include in the output. Supported file formats include VCF files... 
Samples will be outputted in the same order as in the sample file including outputting samples with missing values if no
phenotype information was found in the phenotype files.
Default is to output all samples in the input files in the order they are encountered.
"""


# Notes and TODOs:
# Add option to control what is treated as a missing value
# 	Something like... By default 'NA' is used for missing values in the sample file. Any value in any column that is equal (as a string literal) to "NA" will be treated as missing. (The option -missing-code can be used to alter what is treated as a missing value.) 

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

#import argparse
#import sys
#sys.path.append("/home/fls530/python/mylib")
import click
#import pkclick
import pkcsv as csv
import pkpheno as Pheno
import sys
from pkclick import CSV, SampleList

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  RUN: Commands  %%--

class StdCommand(click.Command):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self.params.insert(0, click.Argument(['files'], nargs=-1, type=click.File()))
		self.epilog = EPILOG_Legal

@click.group()
@click.version_option(version=Version)
def main():
	"""THIS SCRIPT IS STILL EXPERIMENTAL; USE WITH CAUTION"""
	pass

@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to extract.")
def extract(files, columns):
	"""Extract..."""
	pheno = Pheno.Phenotype(csv.DictReader(files[0]), columns=columns)
	for fileobj in files[1:]:
		pheno_new = Pheno.Phenotype(csv.DictReader(fileobj), columns=columns)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()

@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to extract.")
def plink2(files, columns):
	"""Prepare and output phenotupes in psam format for use with Plink2."""
	pheno = Pheno.Psam(csv.DictReader(files[0]), columns=columns)
	for fileobj in files[1:]:
		pheno_new = Pheno.Psam(csv.DictReader(fileobj), columns=columns)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()


@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--covariates', type=CSV(), default="", help="Comma separated list of columns with covariates. Print only these columns (plus any mandatory columns)")
@click.option('-p', '--phenotypes', type=CSV(), default="", help="Comma separated list of columns with phenotypes.")
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTION_samples)
def snptest(files, covariates, phenotypes, samples):
	"""Prepare and output phenotypes in sample file format suitable for use with Snptest."""
	pheno = Pheno.Snptest(csv.DictReader(files[0]), covariates=covariates, phenotypes=phenotypes, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.Snptest(csv.DictReader(fileobj), covariates=covariates, phenotypes=phenotypes, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()

@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to perform normalization on.")
def transform(files, columns):
	"""Used for testing and development."""
	pheno = Pheno.Phenotype(csv.DictReader(files[0]))
	for fileobj in files[1:]:
		pheno_new = Pheno.Phenotype(csv.DictReader(fileobj))
		pheno = pheno.combine_first(pheno_new)
	pheno = pheno.columns_rankINT(columns=columns) # Perform rankINT
	pheno.write()

# --%%  END: Commands  %%--
#
##################################################




##################################################
#
# --%%  RUN: Main Program  %%--

# Execute Main Program

if __name__ == '__main__':
	main()

# --%%  END: Main program  %%--
#
##################################################

