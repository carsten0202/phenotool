

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
from eastwood.eastwood import Incidence, Prevalence
from ukbiobank.ukbiobank import UKBioBank



###########################################################
#
# -%  Click commands for Eastwood algorithms  %-

#
# -%  Incidence Command (Chain Version)  %-

@click.command(name="incidence", no_args_is_help=True)
@click.pass_context
@click.option('-b', '--baseline', default=str(Incidence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-d', '--datediag', default="41270,41280", show_default=True, type=CSV(), help=OPTIONS.datediag)
@click.option('-e', '--enddate', default=str(Incidence.UKBenddate.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.stopdate)
@click.option('-i', '--interval', help=OPTIONS.inciinterval)
@click.option('-p', '--prefix', default="Incidence", show_default=True, help=OPTIONS.columnprefix)
def incidence_ukb(ctx, baseline, prefix, datediag, enddate, interval):
	"""Incidence (diabetes) algorithm from Eastwood2016.
              default=str(date.today()))

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	def processor(pheno):
		incidence = Incidence(pheno, baseline=baseline, datediag=datediag, enddate=enddate, interval=interval)
		pheno[prefix + '_t1dm'] = incidence.t1dm
		pheno[prefix + '_t2dm'] = incidence.t2dm
		pheno[prefix + '_anydm'] = incidence.anydm
		return pheno

	for field in datediag:
		if field not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)
	return processor
	

#
# -%  Prevalence Command (Chain Version)  %-

@click.command(name="prevalence", no_args_is_help=True)
@click.pass_context
@click.option('-a', '--agediag', default="2976,20009", show_default=True, type=CSV(), help=OPTIONS.agediag)
@click.option('-b', '--baseline', default=str(Prevalence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-d', '--datediag', default="41270,41280", show_default=True, type=CSV(), help=OPTIONS.datediag)
@click.option('-e', '--ethnicity', default="21000", show_default=True, type=CSV(), help=OPTIONS.ethnicity)
@click.option('-n', '--name', default="Prevalence", show_default=True, help=OPTIONS.columnname)
@click.option('-r', '--reported', default="2443,4041,20002", show_default=True, type=CSV(), help=OPTIONS.reported)
@click.option('-s', '--style', default="eastwood", show_default=True, type=click.Choice(Prevalence.styles, case_sensitive=False), help=OPTIONS.prevstyles)
@click.option('-t', '--treatments', default="2986,6153,6177,20003", show_default=True, type=CSV(), help=OPTIONS.treatments)
def prevalence_ukb(ctx, agediag, baseline, datediag, ethnicity, name, reported, style, treatments):
	"""Prevalence (diabetes) algorithm from Eastwood2016.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	def processor(pheno):
		prevalence = Prevalence(pheno, agediag=agediag, baseline=baseline, high=datediag, ethnicity=ethnicity, moderate=reported, style=style, treatments=treatments)
		pheno[name] = prevalence.prevalence
		return pheno

	for field in agediag + datediag + ethnicity + reported + treatments:
		if field not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)
	return processor




#
# -%  Incidence Command (Stand-alone Version)  %-

@click.group(chain=True, no_args_is_help=True)
@click.pass_context
@click.option('-b', '--baseline', default=str(Incidence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-d', '--datediag', default="41270,41280", show_default=True, type=CSV(), help=OPTIONS.datediag)
@click.option('-e', '--enddate', default=str(Incidence.UKBenddate.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.stopdate)
@click.option('-i', '--interval', help=OPTIONS.inciinterval)
@click.option('-p', '--prefix', default="Incidence", show_default=True, help=OPTIONS.columnprefix)
def incidence(ctx, baseline, datediag, enddate, interval, prefix):
	"""Incidence (diabetes) algorithm from Eastwood2016.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	# Ensure that ctx.obj exists and is a dict 
	ctx.ensure_object(dict)
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['phenovars'] = datediag
	ctx.obj['to_be_deleted'] = ctx.obj['phenovars']

@incidence.resultcallback()
@click.pass_context
def incidence_pipeline(ctx, processors, baseline, datediag, enddate, interval, prefix):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")

	pheno = ctx.obj['pheno']
	incidence = Incidence(pheno, baseline=baseline, datediag=datediag, enddate=enddate, interval=interval)
	pheno[prefix + '_t1dm'] = incidence.t1dm
	pheno[prefix + '_t2dm'] = incidence.t2dm
	pheno[prefix + '_anydm'] = incidence.anydm

	for processor in processors:
		pheno = processor(pheno)


# CSV output command (Cahined version)
incidence.add_command(Phenotool.csv_chain)

# Plink output command (Chained version)
incidence.add_command(Phenotool.plink_chain)

# CSV output command (Cahined version)
incidence.add_command(Phenotool.rvtest_chain)

# Plink output command (Chained version)
incidence.add_command(Phenotool.snptest_chain)





#
# -%  Prevalence Command (Stand-alone Version)  %-

@click.group(chain=True, invoke_without_command=True, no_args_is_help=True)
@click.pass_context
@click.option('-a', '--agediag', default="2976,20009", show_default=True, type=CSV(), help=OPTIONS.agediag)
@click.option('-b', '--baseline', default=str(Prevalence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-n', '--name', default="Prevalence", show_default=True, help=OPTIONS.columnname)
@click.option('-d', '--datediag', default="41270,41280", show_default=True, type=CSV(), help=OPTIONS.datediag)
@click.option('-e', '--ethnicity', default="21000", show_default=True, type=CSV(), help=OPTIONS.ethnicity)
@click.option('-r', '--reported', default="2443,4041,20002", show_default=True, type=CSV(), help=OPTIONS.reported)
@click.option('-s', '--style', default="eastwood", show_default=True, type=click.Choice(Prevalence.styles, case_sensitive=False), help=OPTIONS.prevstyles)
@click.option('-t', '--treatments', default="2986,6153,6177,20003", show_default=True, type=CSV(), help=OPTIONS.treatments)
def prevalence(ctx, agediag, baseline, datediag, ethnicity, name, reported, style, treatments):
	"""Prevalence (diabetes) algorithm from Eastwood2016

CITATION:
Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	# Ensure that ctx.obj exists and is a dict 
	ctx.ensure_object(dict)
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['phenovars'] = agediag + datediag + ethnicity + reported + treatments
	ctx.obj['to_be_deleted'] = ctx.obj['phenovars']
	ctx.obj['eastwood'] = dict(zip(['agediag', 'baseline', 'ethnicity', 'high', 'moderate', 'style', 'treatments'], [agediag, baseline, ethnicity, datediag, reported, style, treatments]))

@prevalence.resultcallback()
@click.pass_context
def prevalence_pipeline(ctx, processors, agediag, baseline, datediag, ethnicity, name, reported, style, treatments):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")

	pheno = ctx.obj['pheno']
	prevalence = Prevalence(pheno, **ctx.obj.get('eastwood', {}))
	pheno[name] = prevalence.prevalence

	for processor in processors:
		pheno = processor(pheno)


# CSV output command (Cahined version)
prevalence.add_command(Phenotool.csv_chain)

# Plink output command (Chained version)
prevalence.add_command(Phenotool.plink_chain)

# CSV output command (Cahined version)
prevalence.add_command(Phenotool.rvtest_chain)

# Plink output command (Chained version)
prevalence.add_command(Phenotool.snptest_chain)

