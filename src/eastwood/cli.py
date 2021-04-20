

###########################################################
#
# --%%  Setup and Initialize  %%--

__version__ = "0.1"

import click
#import copy
import logging
#import numpy as np
#import pandas as pd
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import phenotool.options as OPTIONS
from pklib.pkclick import CSV, gzFile
import pklib.pkcsv as csv
from eastwood.eastwood import Eastwood
from ukbiobank.ukbiobank import UKBioBank




###########################################################
#
# -%  Click commands for Eastwood algorithms  %-

@click.command()
#@click.pass_context
@click.argument('fileobj', nargs=1, type=gzFile(mode='rb'), metavar="FILE")
@click.option('-a', '--agediag', default="2976,20009", show_default=True, type=CSV(), help=OPTIONS.agediag)
@click.option('-d', '--diagnosed', default="", show_default=True, type=CSV(), help=OPTIONS.diagnosed)
@click.option('-e', '--ethnicity', default="21000", show_default=True, type=CSV(), help=OPTIONS.ethnicity)
@click.option('-r', '--reported', default="2443,4041,20002", show_default=True, type=CSV(), help=OPTIONS.reported)
@click.option('-t', '--treatments', default="2986,6153,6177,20003", show_default=True, type=CSV(), help=OPTIONS.treatments)
def prevalence(fileobj, agediag, diagnosed, ethnicity, reported, treatments):
	"""For Testing do not use."""

# For an independent call we may need this:
#	# Ensure that ctx.obj exists and is a dict 
#	ctx.ensure_object(dict)
	datafields = agediag + diagnosed + ethnicity + reported + treatments
	pheno = UKBioBank(csv.DictReader(fileobj), datafields=datafields)
	eastwood = Eastwood(pheno, agediag=agediag, high=diagnosed, ethnicity=ethnicity, moderate=reported, treatments=treatments)
	pheno['Prevalence'] = eastwood.prevalence
	pheno[['SEX','Prevalence']].write()



#
# -%  Prevalence Command (UKBioBank Version)  %-

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
		if field not in ctx.obj['datafields']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['datafields'].append(field)
	def processor(pheno):
		eastwood = Eastwood(pheno, agediag=agediag, high=diagnosed, ethnicity=ethnicity, moderate=reported, treatments=treatments)
		pheno['Prevalence'] = eastwood.prevalence
		return pheno
	return processor




