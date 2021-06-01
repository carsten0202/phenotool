#!/home/fls530/miniconda3/envs/myscripts/bin/python

###########################################################
#
# ---%%%  Phenotool: Phenotool main file  %%%---
#

# Notes and TODOs:
# Add option to control what is treated as a missing value
# 	Something like... By default 'NA' is used for missing values in the sample file. Any value in any column that is equal (as a string literal) to "NA" will be treated as missing. (The option -missing-code can be used to alter what is treated as a missing value.) 

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

__version__ = """0.12"""

import click
from collections import namedtuple
import logging
import pathlib
import sys

ScriptPath = str(pathlib.Path(__file__).resolve().parent.absolute())
sys.path = [ScriptPath + '/..'] + sys.path

from phenotool.stdcommand import StdCommand
from phenotool.cli import plink_chain, rvtest_chain, snptest_chain, textfile_chain
import phenotool.epilog as EPILOGS
import phenotool.options as OPTIONS
import pklib.pkcsv as csv
from pklib.pkclick import CSV, gzFile, SampleList
from pkpep import pkpep
import ukbiobank.cli as ukbiobank

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  RUN: Commands  %%--

@click.group()
@click.option('--log', default="warning", help=OPTIONS.log, show_default=True)
@click.version_option(version=__version__)
def main(log):
	"""Organize sample information for GWAS analyses.

Read column-based sample information and perform simple sorting, filtering and transformations on phenotype values.
Outputs sample information in formats appropriate for popular GWAS tools including Snptest, RVtest and Plink. 
"""
	try: log_num = getattr(logging, log.upper())
	except AttributeError:
		raise ValueError(f"Invalid log level: '{log}'")
	logging.basicConfig(level=log_num)


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


@main.command(no_args_is_help=True)
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-f', '--fam', type=str, metavar='COLUMN', default=None, help=OPTIONS.fam)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
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
	pheno = Pheno.Psam(csv.DictReader(files[0]), phenovars=columns, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.Psam(csv.DictReader(fileobj), phenovars=columns, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write(header = False if fam else True)


@main.command(no_args_is_help=True)
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
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
	pheno = Pheno.RVtest(csv.DictReader(files[0]), phenovars=columns, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.RVtest(csv.DictReader(fileobj), phenovars=columns, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()


@main.command(no_args_is_help=True)
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--covariates', type=CSV(), default="", help=OPTIONS.covariates)
@click.option('-p', '--phenotypes', type=CSV(), default="", help=OPTIONS.phenotypes)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
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
	pheno = Pheno.Snptest(csv.DictReader(files[0]), covariates=covariates, phenovars=phenotypes, samples=samples)
	for fileobj in files[1:]:
		pheno_new = Pheno.Snptest(csv.DictReader(fileobj), covariates=covariates, phenovars=phenotypes, samples=samples)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()



#
# -%  UKBioBank Command Group  %-

main.add_command(ukbiobank.ukbiobank)



#
# -%  Derive Command Group  %-

@main.group(chain=True)
@click.pass_context
def derive(ctx):
	"""Derive new columns from values in existing ones.

This command group provides access to some of the mathematical functions in the pandas library.
"""
	# ensure that ctx.obj exists and is a dict
	ctx.ensure_object(dict)

@derive.resultcallback()
@click.pass_context
def derive_pipeline(ctx, processors):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")
	pheno = ctx.obj['pheno']
	for processor in processors:
		pheno = processor(pheno)

# Currently implementing transform routines: rankINT is functional
@derive.command(name='rankinv', no_args_is_help=True)
@click.pass_context
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to perform normalization on.")
@click.option('-p', '--prefix', default="rankinv_", show_default=True, help=OPTIONS.columnprefix)
def rankinv_chain(ctx, columns, prefix):
	"""Perform Rank-Based Inverse Normal Transformations."""
	def processor(pheno):
		pheno[[f"{prefix}{col}" for col in columns]] = pheno.columns_rankINT(columns=columns)
		return pheno

	for col in columns:
		if col not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)


# Plink output command (Chained version)
derive.add_command(plink_chain)

# RVtest output command (Chained version)
derive.add_command(rvtest_chain)

# Snptest output command (Chained version)
derive.add_command(snptest_chain)

# Textfile output command (Chained version)
derive.add_command(textfile_chain)

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

