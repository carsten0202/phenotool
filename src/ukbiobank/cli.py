#!/usr/bin/python3

# -*- coding: utf-8 -*-

# TODO:
#               Option to provide a reg ex?
#               Maybe a function to dig up most recent annotation for a given application
#                               ie you'd supply an application nunmber rather than a file?

__version__ = """0.8"""
# v0.6: Fixed some bugs and added textfile as a default command.
# v0.7: Fixed more bugs and added speedup using pd.read_csv.
# v0.8: Added --nrows and fixed more bugfs.

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import eastwood.cli as Eastwood
from phenotool import EPILOG, OPTIONS, plink_chain, rvtest_chain, snptest_chain, textfile_chain
import ukbiobank.options as OPTIONS_UKB
import pklib.pkcsv as csv
from pklib.pkclick import CSV, gzFile, SampleList
from ukbiobank.ukbiobank import UKBioBank

@click.group(chain=True, invoke_without_command=True, no_args_is_help=True, epilog=EPILOG.chained)
@click.option('-d', '--datafields', type=CSV(), required=True, help=OPTIONS_UKB.datafields)
@click.option('-i', '--instances', type=CSV(), help=OPTIONS_UKB.instances)
@click.option('--log', default="warning", show_default=True, help=OPTIONS.log)
@click.option('--nrows', type=int, default=None, help=OPTIONS.nrows)
@click.option('--phenotype-file', type=gzFile(mode='rb'), default=None, envvar='UKBIOBANK_PHENOTYPE_FILE', help=OPTIONS_UKB.phenotype_file)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--sex', type=click.Choice(['genetic','registry'], case_sensitive=False), default='registry', show_default=True, help=OPTIONS_UKB.sex, hidden=True)
@click.option('-v', '--values', default=".", help=OPTIONS.values, hidden=True)
@click.version_option(version=__version__)
@click.pass_context
def ukbiobank(ctx, datafields, instances, log, nrows, phenotype_file, samples, sex, values):
    """UNDER DEVELOPMENT - NOT WELL TESTED.

The UK Biobank lists thousands of phenotypic variables in a strict three tier hierarchy, which requires some
familiarity to navigate effectively. The top of the hierarchy is always the data field which are numerically named.
Most datafields are then divided into instances which typically represent different time points at which data was
colleced. Below the instances we find the data arrays, which are always numerically indexed, but otherwise vary
considerably from datafield to datafield. It is therefore recommended to study the data field descriptions
carefully, for example by using the showcase link given below.

To learn more about the individual datafields their data and their encodings, please refer to:
https://biobank.ndph.ox.ac.uk/showcase/search.cgi
"""
# Would it be relevant to choose between reported and genetic sex?
# 22001-0.0             488264  Categorical (single)        Genetic sex; Uses data-coding 9
# TODO: Fix the hidden options, ie sex & values
# ensure that ctx.obj exists and is a dict
    try: log_num = getattr(logging, log.upper())
    except AttributeError:
        raise ValueError(f"Invalid log level: '{log}'")
    logging.basicConfig(level=log_num)
    ctx.ensure_object(dict)
    if instances:
        datafields = [f"{df}\D{inst}" for df in datafields for inst in instances]
    try: ctx.obj['args']['nrows'] = nrows
    except KeyError:
        ctx.obj['args'] = {'nrows': nrows}
    ctx.obj['args']['phenovars'] = datafields
    ctx.obj['args']['samples'] = samples
    ctx.obj['args']['sexcol'] = sex
    ctx.obj['args']['values'] = values
    ctx.obj['constructor'] = UKBioBank

    # If no subcommand is provided, try the default command
    if ctx.invoked_subcommand is None:
        if phenotype_file is not None:
            result = ctx.invoke(textfile_chain, files=[phenotype_file])
            process_pipeline([result], datafields, instances, log, nrows, phenotype_file, samples, sex, values)
        else:
            logger.error("UKBioBank: You must specify either a command or give a phenotype file. Exiting...")
            exit(1)

@ukbiobank.result_callback()
@click.pass_obj
def process_pipeline(obj, processors, datafields, instances, log, nrows, phenotype_file, samples, sex, values):
    logger.debug(f"Pipeline: Cols to be deleted: {obj.get('to_be_deleted')}")
    logger.debug(f'Processors: {processors}')
    pheno = obj['pheno']
    for processor in processors:
        pheno = processor(pheno)



#
# -%    Add Command on External Commands (Chained Versions)  %-

# CSV output command (Chained version)
ukbiobank.add_command(textfile_chain)

# Eastwood Incidence Command (UKBioBank Version)
ukbiobank.add_command(Eastwood.incidence_ukb)

# Plink output command (Chained version)
ukbiobank.add_command(plink_chain)

# Eastwood Prevalence Command (UKBioBank Version)
ukbiobank.add_command(Eastwood.prevalence_ukb)

# RVtest output command (Chained version)
ukbiobank.add_command(rvtest_chain)

# Snptest output command (Chained version)
# ukbiobank.add_command(snptest_chain)
# TODO: SNPTEST IS CURRENTLY BROKEN

