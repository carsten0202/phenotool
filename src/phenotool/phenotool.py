#!/home/fls530/miniconda3/envs/myscripts/bin/python

# phenotool.py

# Notes and TODOs:
# Add option to control what is treated as a missing value
# 	Something like... By default 'NA' is used for missing values in the sample file. Any value in any column that is equal (as a string literal) to "NA" will be treated as missing. (The option -missing-code can be used to alter what is treated as a missing value.) 

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

__version__ = """0.8 (Development Version)"""

import click
from collections import namedtuple
import logging
import pathlib
import sys

ScriptPath = str(pathlib.Path(__file__).resolve().parent.absolute())
sys.path = [ScriptPath + '/..'] + sys.path
sys.path = [ScriptPath + '/../src'] + sys.path

import tmp.pkcsv as csv
from pkclick import CSV, SampleList


EPILOG = namedtuple('Epilog', ['legal'])(
legal = """
Written by Carsten Friis Rundsten <carsten.rundsten@sund.ku.dk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
""")

OPTION = namedtuple('Options', ['columns','covariates','fam','files','log','phenotypes','samples'])(
columns = """
Comma separated list of columns to output in addition to mandatory columns. Default is to output all columns.
""",
covariates = """
Comma separated list of columns with covariates. Print only these columns (plus any mandatory columns)
""",
fam = """
Output fam format for plink1.9. If specified, takes one mandatory argument stating the column name of the phenotype to
include in the sixth fam file column.
""",
files = """
Input File(s). %(prog)s accepts one or more input files including '-' symbolizing stdin. The precise format of each
input file will be autodetected, but it should be some form of delimited text data file like 'csv' or tab-delimited.
""",
log = """
Control logging. Valid levels: 'debug', 'info', 'warning', 'error', 'critical'.
""",
phenotypes = """
Comma separated list of columns with phenotypes.
""",
samples="""
File with samples to include in the output. Samples will be outputted in the exact same order as in the sample file
including outputting samples with missing values if no phenotype information was found. The sample file can be a plain text file with sample names or a VCF file with sample genotypes.
""",
)

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
		self.epilog = EPILOG.legal

@click.group()
@click.version_option(version=__version__)
@click.option('--log', default="warning", help=OPTION.log, show_default=True)
def main(log):
	"""Organize sample information for GWAS analyses.

Read column-based sample information and perform simple sorting, filtering and transformations on phenotype values.
Outputs sample information in formats appropriate for popular GWAS tools including Snptest, RVtest and Plink. 
"""
	try:
		log_num = getattr(logging, log.upper())
	except AttributeError:
		raise ValueError(f"Invalid log level: '{log}'")
	logging.basicConfig(level=log_num)


#@main.command(cls=StdCommand, no_args_is_help=True)
#@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to extract.")
#def extract(files, columns):
#	"""Extract not implemented yet. Don't really know if we need this guy, or if the options to extract cannot be handled thorugh the other commands?"""
#	pheno = Pheno.Phenotype(csv.DictReader(files[0]), columns=columns)
#	for fileobj in files[1:]:
#		pheno_new = Pheno.Phenotype(csv.DictReader(fileobj), columns=columns)
#		pheno = pheno.combine_first(pheno_new)
#	pheno.write()

from pkpep import pkpep

main.add_command(pkpep.main, name="pep")
#@main.command(cls=StdCommand, no_args_is_help=True, hidden=True)
#@click.option('-c', '--columns', type=CSV(), default="", help=OPTION.columns)
#def pep(files, columns):
"""NOT IMPLEMENTED YET. Output CSV file suitable for PEP.

The Portable Encapsulated Projects (PEP for short) community effort to facilitate the portability, reusability and
durability of sample metadata.

\b
For more on the PEP community effort, please refer to:
http://pep.databio.org/en/latest/
"""
#	import pkpheno as Pheno
#	assert True, "."
#	pheno = Pheno.PEP(csv.DictReader(files[0]), columns=columns)
#	for fileobj in files[1:]:
#		pheno_new = Pheno.PEP(csv.DictReader(fileobj), columns=columns)
#	pheno = pheno.combine_first(pheno_new)
#	pheno.write()


@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--columns', type=CSV(), default="", help=OPTION.columns)
@click.option('-f', '--fam', type=str, metavar='COLUMN', default=None, help=OPTION.fam)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTION.samples)
def plink(files, columns, fam, samples):
	"""Output phenotypes in psam/fam format for use with Plink.

A properly formatted psam file has in addition to the phenotype columns one or more of the following recognizable
columns: 'IID' (individual ID; required), 'FID', 'SID', 'PAT', 'MAT' and 'SEX'. If no columns with those exact names
are present in the input, then the program tries to guess which input columns might contain the missing information.
The guessing uses common synonyms, eg mapping 'gender' to 'SEX'. It also maps across supported formats mapping 'ID_1'
(snptest) to 'IID'.

For more on the psam format please refer to:
https://www.cog-genomics.org/plink/2.0/formats#psam

A properly formatted fam file (plink1.9) has *no* header, but expects the following six columns in exact order: 'FID',
'IID', 'PAT', MAT', 'SEX', and one final phenotype column. The plink1.9 fam format only supports one phenotype. To work
with more than one phenotype in plink1.9 the psam files prepared by this program are designed to be readable as an
alternate phenotype file in plink1.9.

\b
For more on alternate phenotype files in plink1.9, please refer to:
https://www.cog-genomics.org/plink/1.9/input#pheno

\b
For more on the fam format please refer to:
https://www.cog-genomics.org/plink/1.9/formats#fam
"""
	import pkpheno as Pheno
	assert sum([1 for x in [columns,fam] if x]) <= 1, "'--columns' and '--fam' are mutually exclusive; please only specify one of them."
	if fam:
		columns = [fam]
	pheno = Pheno.Psam(csv.DictReader(files[0]), columns=columns, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.Psam(csv.DictReader(fileobj), columns=columns, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write(header = False if fam else True)


@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--columns', type=CSV(), default="", help=OPTION.columns)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTION.samples)
def rvtest(files, columns, samples):
	"""UNTESTED; Output phenotypes in psam-like format for RVtest.

RVtest phenotype files are very similar to the psam format. They are essentially plink2 files with a few caveats. The
names 'fatid' and 'matid' are used for paternal and maternal ids and 'sex' is encoded 0=males,1=females.

\b
For more on RVtest phenotype files, please refer to:
http://zhanxw.github.io/rvtests/#phenotype-file
"""
	import pkpheno as Pheno
	import pandas as pd
	pheno = Pheno.RVtest(csv.DictReader(files[0]), columns=columns, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.RVtest(csv.DictReader(fileobj), columns=columns, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()


@main.command(cls=StdCommand, no_args_is_help=True)
@click.option('-c', '--covariates', type=CSV(), default="", help="Comma separated list of columns with covariates. Print only these columns (plus any mandatory columns)")
@click.option('-p', '--phenotypes', type=CSV(), default="", help="Comma separated list of columns with phenotypes.")
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTION.samples)
def snptest(files, covariates, phenotypes, samples):
	"""Output phenotypes in sample format for use with Snptest.

A properly formatted snptest sample (*.sam) file must contain the columns 'ID_1', 'ID_2', 'missing' and 'sex' in that
order, followed by any covariate columns and finally by any phenotype columns. If no columns with those exact names are
present in the input, then the program first tries to guess which input columns might contain the missing information.
The guessing uses common synonyms, eg mapping 'gender' to 'sex'. It also maps across supported formats mapping 'IID'
(Plink) to 'ID_1'. If the guessing fails, it then tries to fill in the missing information, which will produce
functional files for columns 'ID_2' (using 'ID_1' values) and 'missing' (using missing values).

Furthermore a snptest sample file must have phenotypes and covariates explicitly stated as such. Since this information
cannot be inferred form the data alone, the user should provide this using the '--covariates' and '--phenotypes'
options documented below.

\b
For the official docs on the sample format please refer to:
https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html

\b
Unofficial, but good (Scroll down):
https://jmarchini.org/file-formats/
"""
	import pkpheno as Pheno
	pheno = Pheno.Snptest(csv.DictReader(files[0]), covariates=covariates, phenotypes=phenotypes, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.Snptest(csv.DictReader(fileobj), covariates=covariates, phenotypes=phenotypes, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()

# Currently implementing transform routines: rankINT is functional
@main.command(cls=StdCommand, no_args_is_help=True, hidden=True)
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to perform normalization on.")
def test(files, columns):
	"""Used for testing and development."""
	import pkpheno as Pheno
	pheno = Pheno.Phenotype(csv.DictReader(files[0]))
	for fileobj in files[1:]:
		pheno_new = Pheno.Phenotype(csv.DictReader(fileobj))
		pheno = pheno.combine_first(pheno_new)
	pheno = pheno.columns_rankINT(columns=columns) # Perform rankINT
	pheno.write()

# --%%  END: Commands  %%--
#
##################################################



from ukbiobank import ukbiobank

main.add_command(ukbiobank.ukbiobank)



##################################################
#
# --%%  RUN: Main Program  %%--

# Execute Main Program

if __name__ == '__main__':
	main()

# --%%  END: Main program  %%--
#
##################################################

