
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

__version__ = """0.12.1"""
# v0.12.1: Fixed a ton of parsing errors

import click
from collections import namedtuple
import logging
import pathlib
import sys

ScriptPath = str(pathlib.Path(__file__).resolve().parent.absolute())
sys.path = [ScriptPath + '/..'] + sys.path

from phenotool.stdcommand import StdCommand
from phenotool import OPTIONS, plink, plink_chain, rvtest, rvtest_chain, snptest, snptest_chain, textfile, textfile_chain
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
@click.pass_context
@click.option('--log', default="warning", help=OPTIONS.log, show_default=True)
@click.version_option(version=__version__)
def main(ctx, log):
    """Organize sample information for GWAS analyses.

Read column-based sample information and perform simple sorting, filtering and transformations on phenotype values.
Outputs sample information in formats appropriate for popular GWAS tools including Snptest, RVtest and Plink. 
"""
    try: log_num = getattr(logging, log.upper())
    except AttributeError:
        raise ValueError(f"Invalid log level: '{log}'")
    logging.basicConfig(level=log_num)
    ctx.ensure_object(dict)


# PEP Command
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

# Plink Command
main.add_command(plink)

# RVtest Command
main.add_command(rvtest)

# Snptest Command
#main.add_command(snptest)
# TODO: Fix Snptest (Broken with set_columns)

# TextFile Command
main.add_command(textfile)

# UKBioBank Command Group
main.add_command(ukbiobank.ukbiobank)



#
# -%  Derive Command Group  %-

@main.group(chain=True)
@click.pass_obj
def derive(obj):
    """Derive new columns from values in existing ones.

This command group provides access to some of the mathematical functions in the pandas library.
"""
    try: obj['args']
    except KeyError: obj['args'] = dict()
    try: obj['args']['phenovars']
    except KeyError: obj['args']['phenovars'] = list()

@derive.result_callback()
@click.pass_obj
def derive_pipeline(obj, processors):
    logging.debug(f"Pipeline: Cols to be deleted: {obj.get('to_be_deleted')}")
    pheno = obj['pheno']
    for processor in processors:
        pheno = processor(pheno)

# Currently implementing transform routines: rankINT is functional
@derive.command(name='rankinv', no_args_is_help=True)
@click.pass_obj
@click.option('-c', '--columns', type=CSV(), default="", help="Comma separated list of columns to perform normalization on.")
@click.option('-p', '--prefix', default="rankinv_", show_default=True, help=OPTIONS.columnprefix)
def rankinv_chain(obj, columns, prefix):
    """Perform Rank-Based Inverse Normal Transformations."""
    def processor(pheno):
        pheno[[f"{prefix}{col}" for col in columns]] = pheno.derive_rankINT(columns=columns)
        return pheno

    for col in columns:
        if col not in obj['args']['phenovars']:
            obj['to_be_deleted'] = obj.get('to_be_deleted', list()) + [col]
            obj['args']['phenovars'].append(col)
    return processor


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

