#!/usr/bin/python3

# -*- coding: utf-8 -*-

# TODO:
#       Option to provide a reg ex?
#       Maybe a function to dig up most recent annotation for a given application
#               ie you'd supply an application nunmber rather than a file?

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import eastwood.cli as eastwood
import phenotool.options as OPTIONS
import pklib.pkcsv as csv
from pklib.pkclick import CSV, gzFile, SampleList
from pkpheno.pkpheno import Phenotype
from ukbiobank.ukbiobank import UKBioBank

@click.group(chain=True, invoke_without_command=True, no_args_is_help=True)
# @click.argument('fileobj', nargs=1, type=gzFile(mode='rb'), metavar="FILE")
@click.option('-d', '--datafields', required=True, type=CSV(), help="Data Field(s) to output. Several fields can be specified as a comma-separated string with no spaces.")
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('-v', '--value', default=".", help="Only subject having this value in at least one of the specified data field(s) will be printed. Note that this uses text-based Regular Expression matching, so '-v 1' will match a value of '12'. Use '-v \"\\b1\\b\"' to match only '1'. Default is to print all subjects having at least one non-empty value in the specified data field(s).")
@click.pass_context
def ukbiobank(ctx, datafields, samples, value):
	"""UNDER DEVELOPMENT - DO NOT USE."""
	logger.debug(f"From samplefile: {samples}")

# Would it be relevant to choose between reported and genetic sex?
# 22001-0.0       488264  Categorical (single)    Genetic sex; Uses data-coding 9

	# ensure that ctx.obj exists and is a dict
	ctx.ensure_object(dict)
	ctx.obj['datafields'] = datafields
	ctx.obj['samples'] = samples
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['DEBUG'] = True



@ukbiobank.resultcallback()
@click.pass_context
def process_pipeline(ctx, processors, datafields, samples, value):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")
	pheno = ctx.obj['pheno']
	for processor in processors:
		pheno = processor(pheno)



#
# -%  Eastwood Prevalence Command (UKBioBank Version)  %-

ukbiobank.add_command(eastwood.prevalence_ukb)


#
# -%  Output Commands (Chained)  %-

@ukbiobank.command(name="csv")
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
def csvfile(ctx, files):
	"""Output as CSV."""
	ctx.obj['constructor'] = ctx.obj.get('constructor', Phenotype)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), datafields=ctx.obj['datafields'], samples=ctx.obj['samples'])
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), datafields=ctx.obj['datafields'], samples=ctx.obj['samples'])
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno.write()
	return processor




