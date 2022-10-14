
###########################################################
#
# ---%%%  Phenotool: cli.py Command Line Interface  %%%---
#

#
# -%  Setup  %-

import click

import phenotool.options as OPTIONS
from pklib.pkclick import CSV, gzFile, SampleList
import pklib.pkcsv as csv



###########################################################
#
# --%%  Chained Output Commands (Those which must be last)  %%--

#
# -%  Plink Command; Chained Version  %-

@click.command(name="plink", no_args_is_help=True)
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-f', '--fam', type=str, metavar='COLUMN', default=None, help=OPTIONS.fam)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def plink_chain(ctx, files, columns, fam, samples):
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
	from pkpheno.pkpheno import Psam
	assert sum([1 for x in [columns,fam] if x]) <= 1, "'--columns' and '--fam' are mutually exclusive; please only specify one of them."
	if fam:
		columns = [fam]
		ctx.obj['fam'] = True
	ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + columns)) # Clever little trick to get unique list
	if samples:
		ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples)) if ctx.obj.get('samples') else samples
	ctx.obj['constructor'] = ctx.obj.get('constructor', Psam)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)

	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno = pheno.to_psam()
		pheno.write(header = False if ctx.obj.get('fam') else True)
		return pheno
	return processor



#
# -%  RVtest Command; Chained Version  %-

@click.command(name="rvtest", no_args_is_help=True)
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def rvtest_chain(ctx, files, columns, samples):
	"""UNTESTED; Output phenotypes in psam-like format for RVtest.

RVtest phenotype files are very similar to the psam format. They are essentially plink2 files with a few caveats. The
names 'fatid' and 'matid' are used for paternal and maternal ids and 'sex' is encoded 0=males,1=females.

\b
For more on RVtest phenotype files, please refer to:
http://zhanxw.github.io/rvtests/#phenotype-file
"""
	from pkpheno.pkpheno import RVtest
	ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + columns)) # Clever little trick to get unique list
	if samples:
		ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
	ctx.obj['constructor'] = ctx.obj.get('constructor', RVtest)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)

	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno = pheno.to_rvtest()
		pheno.write()
		return pheno
	return processor



#
# -%  Snptest Command; Chained Version  %-

@click.command(name="snptest", no_args_is_help=True)
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--covariates', type=CSV(), default="", help=OPTIONS.covariates)
@click.option('-p', '--phenotypes', type=CSV(), default="", help=OPTIONS.phenotypes)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def snptest_chain(ctx, files, covariates, phenotypes, samples):
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
	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno = pheno.to_snptest(covariates = covariates)
		pheno.write()
		return pheno

	from pkpheno.pkpheno import Snptest
	ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + covariates + phenotypes)) # Clever little trick to get unique list
	if samples:
		ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
	ctx.obj['constructor'] = ctx.obj.get('constructor', Snptest)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
	return processor



#
# -%  Plain Text File Command; Chained Version  %-

@click.command(name="textfile", no_args_is_help=True)
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--csv', 'formatflag', flag_value='csv', default=True, help=OPTIONS.csv)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--tsv', 'formatflag', flag_value='tsv', help=OPTIONS.tsv)
def textfile_chain(ctx, files, columns, formatflag, samples):
	"""Output phenotypes in customizable text format."""
	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno = pheno.to_textfile()
		pheno.write()
		return pheno

	assert formatflag == 'csv', "Sorry, textfile currently only supports csv output. This will change in future versions."

	ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + columns)) # Clever little trick to get unique list
	if samples:
		ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
	from pkpheno import TextFile
	ctx.obj['constructor'] = ctx.obj.get('constructor', TextFile)
	ctx.obj['pheno'] = ctx.obj['constructor'](files[0], phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](fileobj, phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
	return processor

