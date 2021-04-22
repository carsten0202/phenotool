

import click

import phenotool.options as OPTIONS
from pklib.pkclick import CSV, gzFile, SampleList
import pklib.pkcsv as csv
from pkpheno.pkpheno import Phenotype, Psam



###########################################################
#
# --%%  Chained Output Commands (Those which must be last)  %%--

#
# -%  Plain CSV Command; Chained Version  %-

@click.command(name="csv")
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
def csv_chain(ctx, files):
	"""Output as CSV."""
	ctx.obj['constructor'] = ctx.obj.get('constructor', Phenotype)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj['samples'])
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj['samples'])
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno.write()
	return processor




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
	assert sum([1 for x in [columns,fam] if x]) <= 1, "'--columns' and '--fam' are mutually exclusive; please only specify one of them."
	if fam:
		columns = [fam]
		ctx.obj['fam'] = True
	ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + columns)) # Clever little trick to get unique list
	if samples:
		ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
	ctx.obj['constructor'] = ctx.obj.get('constructor', Psam)
	ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj['samples'])
	for fileobj in files[1:]:
		pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj['samples'])
		ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)

	def processor(pheno):
		if ctx.obj.get('to_be_deleted'):
			pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
		pheno = pheno.to_psam()
		pheno.write(header = False if ctx.obj.get('fam') else True)
		return pheno
	return processor

