#!/usr/bin/python3

# -*- coding: utf-8 -*-

# TODO:
#       Option to provide a reg ex?
#       Maybe a function to dig up most recent annotation for a given application
#               ie you'd supply an application nunmber rather than a file?

__version__ = """0.3.0"""

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import eastwood.cli as Eastwood
import phenotool.cli as Phenotool
import phenotool.options as OPTIONS
import pklib.pkcsv as csv
from pklib.pkclick import CSV, gzFile, SampleList
from pkpheno.pkpheno import Phenotype
from ukbiobank.ukbiobank import UKBioBank

@click.group(chain=True, invoke_without_command=True, no_args_is_help=True)
@click.option('-d', '--datafields', type=CSV(), default="", help=OPTIONS.datafields)
@click.option('-i', '--instances', type=CSV(), help=OPTIONS.instances)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('-v', '--values', default=".", help=OPTIONS.values)
@click.version_option(version=__version__)
@click.pass_context
def ukbiobank(ctx, datafields, instances, samples, values):
	"""UNDER DEVELOPMENT - NOT WELL TESTED.

The UK Biobank lists thousands of phenotypic variables in a strict three tier hierarchy, which requires some
familiarity to navigate effectively. The top of the hierarchy is always the data field which are numerically named.
Most datafields are then divided into instances which typically represent different time points at which data was
colleced. Below the instances we find the data arrays, which are always numerically indexed, but otherwise vary
considerably from datafield to datafield. It is therefore recommended to study the data field descriptions
carefully, for example by using the showcase link given below.

To learn more about the individual datafields, their data and their encodings please refer to:
https://biobank.ndph.ox.ac.uk/showcase/search.cgi
"""
#	logger.debug(f"From samplefile: {samples}")

# Would it be relevant to choose between reported and genetic sex?
# 22001-0.0       488264  Categorical (single)    Genetic sex; Uses data-coding 9

	# ensure that ctx.obj exists and is a dict
	ctx.ensure_object(dict)

	if instances:
#		sys.exit("Sorry, Instances is not implemented yet.")
		datafields = [f"{df}\D{inst}" for df in datafields for inst in instances]

	ctx.obj['phenovars'] = datafields
	ctx.obj['samples'] = samples
	ctx.obj['values'] = values
	ctx.obj['constructor'] = UKBioBank


@ukbiobank.resultcallback()
@click.pass_context
def process_pipeline(ctx, processors, datafields, instances, samples, values):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")
	pheno = ctx.obj['pheno']
	for processor in processors:
		pheno = processor(pheno)



#
# -%  Add Command on External Commands (Chained Versions)  %-

# CSV output command (Chained version)
ukbiobank.add_command(Phenotool.csv_chain)

# Eastwood Incidence Command (UKBioBank Version)
ukbiobank.add_command(Eastwood.incidence_ukb)

# Plink output command (Chained version)
ukbiobank.add_command(Phenotool.plink_chain)

# Eastwood Prevalence Command (UKBioBank Version)
ukbiobank.add_command(Eastwood.prevalence_ukb)

# RVtest output command (Chained version)
ukbiobank.add_command(Phenotool.rvtest_chain)

# Snptest output command (Chained version)
ukbiobank.add_command(Phenotool.snptest_chain)

