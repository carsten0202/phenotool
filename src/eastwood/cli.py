

###########################################################
#
# --%%  Setup and Initialize  %%--

__version__ = "0.4"

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

from phenotool import EPILOG, OPTIONS, plink_chain, rvtest_chain, snptest_chain, textfile_chain
from pklib.pkclick import CSV, gzFile, Timedelta
import pklib.pkcsv as csv
from eastwood.eastwood import Incidence, Prevalence
from ukbiobank.ukbiobank import UKBioBank



###########################################################
#
# -%  Click commands for Eastwood algorithms  %-

#
# -%  Incidence Command (Chain Version)  %-

# default=str(date.today()))
@click.command(name="incidence", no_args_is_help=True, epilog=EPILOG.chained)
@click.pass_context
@click.option('-b', '--baseline', default=str(Incidence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-e', '--enddate', default=str(Incidence.UKBenddate.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.stopdate)
@click.option('-i', '--interval', type=Timedelta(), help=OPTIONS.inciinterval)
@click.option('-p', '--prefix', default="Incidence", show_default=True, help=OPTIONS.columnprefix)
def incidence_ukb(ctx, baseline, prefix, enddate, interval):
	"""Incidence (diabetes) algorithm from Eastwood2016.

Incidence is calculated by first consulting a prevalence algorithm (same as for the 'prevalence' command) to mask out
any cases occuring prior to the baseline date. It then draws upon secondary care data from the UKBiobank, originally
drawn from the English Clinical Practice Research, to calculate incidence rates between the baseline and end dates.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	def processor(pheno):
		incidence = Incidence(pheno, baseline=baseline, enddate=enddate, interval=interval)
		pheno[prefix + '_t1dm'] = incidence.t1dm
		pheno[prefix + '_t2dm'] = incidence.t2dm
		pheno[prefix + '_anydm'] = incidence.anydm
		return pheno

	for field in Incidence.UKBioFields:
		if field not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)
	return processor
	


#
# -%  Prevalence Command (Chain Version)  %-

@click.command(name="prevalence", no_args_is_help=True, epilog=EPILOG.chained)
@click.pass_context
@click.option('-b', '--baseline', default=str(Prevalence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('--date/--no-date', default=False, show_default=True, help="NOT IMPLEMENTED")
@click.option('--datename', default="Prevalence_date", show_default=True, help="NOT IMPLEMENTED")
@click.option('-n', '--name', default="Prevalence", show_default=True, help=OPTIONS.columnname)
@click.option('-s', '--style', default="eastwood", show_default=True, type=click.Choice(Prevalence.styles, case_sensitive=False), help=OPTIONS.prevstyles)
def prevalence_ukb(ctx, baseline, date, datename, name, style):
	"""Prevalence (diabetes) algorithm from Eastwood2016.

The algorithm uses UK Biobank self-reported medical history and medication as well as hospital in-patient data to
assign prevalent diabetes and type.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	def processor(pheno):
		prevalence = Prevalence(pheno, baseline=baseline, style=style)
		pheno[name] = prevalence.prevalence
		if date:
			logger.debug(f"date = {datename}")
			pheno[datename] = prevalence.datediag()
		return pheno

	for field in Prevalence.UKBioFields:
		if field not in ctx.obj['phenovars']:
			ctx.obj['to_be_deleted'] = ctx.obj.get('to_be_deleted', list()) + [field]
			ctx.obj['phenovars'].append(field)
	return processor




#
# -%  Incidence Command (Stand-alone Version)  %-

@click.group(chain=True, no_args_is_help=True)
@click.pass_context
@click.option('-b', '--baseline', default=str(Incidence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-e', '--enddate', default=str(Incidence.UKBenddate.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.stopdate)
@click.option('-i', '--interval', type=Timedelta(), help=OPTIONS.inciinterval)
@click.option('-p', '--prefix', default="Incidence", show_default=True, help=OPTIONS.columnprefix)
@click.version_option(version=__version__)
def incidence(ctx, baseline, enddate, interval, prefix):
	"""Incidence (diabetes) algorithm from Eastwood2016.

Incidence is calculated by first consulting a prevalence algorithm (same as for the 'prevalence' command) to mask out
any cases occuring prior to the baseline date. It then draws upon secondary care data from the UKBiobank, originally
drawn from the English Clinical Practice Research, to calculate incidence rates between the baseline and end dates.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	# Ensure that ctx.obj exists and is a dict 
	ctx.ensure_object(dict)
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['phenovars'] = Incidence.UKBioFields
	ctx.obj['to_be_deleted'] = Incidence.UKBioFields

@incidence.result_callback()
@click.pass_context
def incidence_pipeline(ctx, processors, baseline, enddate, interval, prefix):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")

	pheno = ctx.obj['pheno']
	incidence = Incidence(pheno, baseline=baseline, enddate=enddate, interval=interval)
	pheno[prefix + '_t1dm'] = incidence.t1dm
	pheno[prefix + '_t2dm'] = incidence.t2dm
	pheno[prefix + '_anydm'] = incidence.anydm

	for processor in processors:
		pheno = processor(pheno)


# CSV output command (Cahined version)
incidence.add_command(textfile_chain)

# Plink output command (Chained version)
incidence.add_command(plink_chain)

# CSV output command (Cahined version)
incidence.add_command(rvtest_chain)

# Plink output command (Chained version)
incidence.add_command(snptest_chain)





#
# -%  Prevalence Command (Stand-alone Version)  %-

@click.group(chain=True, invoke_without_command=True, no_args_is_help=True)
@click.pass_context
@click.option('-b', '--baseline', default=str(Prevalence.UKBbaseline.date()), show_default=True, type=click.DateTime(formats=["%Y-%m-%d"]), help=OPTIONS.baseline)
@click.option('-n', '--name', default="Prevalence", show_default=True, help=OPTIONS.columnname)
@click.option('-s', '--style', default="eastwood", show_default=True, type=click.Choice(Prevalence.styles, case_sensitive=False), help=OPTIONS.prevstyles)
@click.version_option(version=__version__)
def prevalence(ctx, baseline, name, style):
	"""Prevalence (diabetes) algorithm from Eastwood2016

The algorithm uses UK Biobank self-reported medical history and medication as well as hospital in-patient data to
assign prevalent diabetes and type.

CITATION:

Algorithms for the capture and adjudication of prevalent and incident diabetes in UK Biobank
SV Eastwood, R Mathur, M Atkinson, S Brophy.
PloS one 2016 Sep 15; 11(9): e0162388
https://doi.org/10.1371/journal.pone.0162388
"""
	# Ensure that ctx.obj exists and is a dict 
	ctx.ensure_object(dict)
	ctx.obj['constructor'] = UKBioBank
	ctx.obj['phenovars'] = Prevalence.UKBioFields
	ctx.obj['to_be_deleted'] = Prevalence.UKBioFields

@prevalence.result_callback()
@click.pass_context
def prevalence_pipeline(ctx, processors, baseline, name, style):
	logger.debug(f"Pipeline: Cols to be deleted: {ctx.obj.get('to_be_deleted')}")

	pheno = ctx.obj['pheno']
	prevalence = Prevalence(pheno, baseline=baseline, style=style)
	pheno[name] = prevalence.prevalence

	for processor in processors:
		pheno = processor(pheno)


# CSV output command (Cahined version)
prevalence.add_command(textfile_chain)

# Plink output command (Chained version)
prevalence.add_command(plink_chain)

# CSV output command (Cahined version)
prevalence.add_command(rvtest_chain)

# Plink output command (Chained version)
prevalence.add_command(snptest_chain)

