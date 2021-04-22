

###########################################################
#
# --%%  Setup and Initialize  %%--

__version__ = "0.2"

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import phenotool.cli as Phenotool
import phenotool.options as OPTIONS
from pklib.pkclick import CSV, gzFile
import pklib.pkcsv as csv
from eastwood.eastwood import Eastwood
from ukbiobank.ukbiobank import UKBioBank



###########################################################
#
# -%  Click commands for Eastwood algorithms  %-

#
# -%  Prevalence Command (Chain Version)  %-

@click.command(name="prevalence")
@click.pass_context
@click.option('-a', '--agediag', default="2976,20009", show_default=True, type=CSV(), help=OPTIONS.agediag)
@click.option('-d', '--diagnosed', default="", show_default=True, type=CSV(), help=OPTIONS.diagnosed)
@click.option('-e', '--ethnicity', default="21000", show_default=True, type=CSV(), help=OPTIONS.ethnicity)
@click.option('-r', '--reported', default="2443,4041,20002", show_default=True, type=CSV(), help=OPTIONS.reported)
@click.option('-t', '--treatments', default="2986,6153,6177,20003", show_default=True, type=CSV(), help=OPTIONS.treatments)
def prevalence_ukb(ctx, agediag, diagnosed, ethnicity, reported, treatments):
	"""The Prevalence algorithm from Eastwood2016"""
	for field in agediag + diagnosed + ethnicity + reported + treatments:
		if field not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)
	def processor(pheno):
		eastwood = Eastwood(pheno, agediag=agediag, high=diagnosed, ethnicity=ethnicity, moderate=reported, treatments=treatments)
		pheno['Prevalence'] = eastwood.prevalence
		return pheno
	return processor



#
# -%  Prevalence Command (Stand-alone Version)  %-

#@click.command(no_args_is_help=True))
@click.group(chain=True, invoke_without_command=True, no_args_is_help=True)
@click.pass_context
@click.option('-a', '--agediag', default="2976,20009", show_default=True, type=CSV(), help=OPTIONS.agediag)
@click.option('-d', '--diagnosed', default="", show_default=True, type=CSV(), help=OPTIONS.diagnosed)
@click.option('-e', '--ethnicity', default="21000", show_default=True, type=CSV(), help=OPTIONS.ethnicity)
@click.option('-r', '--reported', default="2443,4041,20002", show_default=True, type=CSV(), help=OPTIONS.reported)
@click.option('-t', '--treatments', default="2986,6153,6177,20003", show_default=True, type=CSV(), help=OPTIONS.treatments)
def prevalence(ctx, agediag, diagnosed, ethnicity, reported, treatments):
	"""For Testing do not use."""
	# Ensure that ctx.obj exists and is a dict 
	ctx.ensure_object(dict)
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['phenovars'] = agediag + diagnosed + ethnicity + reported + treatments
	ctx.obj['samples'] = None
	ctx.obj['to_be_deleted'] = ctx.obj['phenovars']
	ctx.obj['eastwood'] = dict(zip(['agediag', 'ethnicity', 'high', 'moderate', 'treatments'], [agediag, ethnicity, diagnosed, reported, treatments]))


@prevalence.resultcallback()
@click.pass_context
def process_pipeline(ctx, processors, agediag, diagnosed, ethnicity, reported, treatments):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")

	pheno = ctx.obj['pheno']
	eastwood = Eastwood(pheno, **ctx.obj.get('eastwood', {}))
	pheno['Prevalence'] = eastwood.prevalence

	for processor in processors:
		pheno = processor(pheno)



#
# -%  Chained-Subcommands  %-

# CSV output command (Cahined version)
prevalence.add_command(Phenotool.csv_chain)

# Plink output command (Chained version)
prevalence.add_command(Phenotool.plink_chain)

