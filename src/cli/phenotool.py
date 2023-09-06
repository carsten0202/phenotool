
###########################################################
#
# ---%%%  Phenotool: Phenotool main file  %%%---
#

# Notes and TODO:
# Add option to control what is treated as a missing value
# 	Something like... By default 'NA' is used for missing values in the sample file. Any value in any column that is equal (as a string literal) to "NA" will be treated as missing. (The option -missing-code can be used to alter what is treated as a missing value.) 

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

import click
from collections import namedtuple
import logging
import pathlib
import sys

ScriptPath = str(pathlib.Path(__file__).resolve().parent.absolute())
sys.path = [ScriptPath + '/..'] + sys.path

from cli.version import __version__
from cli.derive import derive
from cli.ukbiobank import ukbiobank

from phenotool.stdcommand import StdCommand
from phenotool import OPTIONS, plink, rvtest, snptest, textfile
# import pklib.pkcsv as csv
# from pklib.pkclick import CSV, isalFile, SampleList
from pkpep import pkpep

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
    """Organize column-based information for analyses.

Read column-based sample information and perform simple sorting, filtering and transformations on phenotype values.
Outputs sample information in formats appropriate for popular GWAS tools including Snptest, RVtest and Plink. 
"""
    try: log_num = getattr(logging, log.upper())
    except AttributeError:
        raise ValueError(f"Invalid log level: '{log}'")
    logging.basicConfig(level=log_num)
    ctx.ensure_object(dict)
    ctx.obj['files'] = ctx.obj.get('files', [])


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

# Derive Command Group 
main.add_command(derive)

# Plink Command
main.add_command(plink)

# RVtest Command
main.add_command(rvtest)

# Snptest Command
main.add_command(snptest)

# TextFile Command
main.add_command(textfile)

# UKBioBank Command Group
main.add_command(ukbiobank)

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

