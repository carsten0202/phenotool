#!/usr/bin/python

import click

from pklib.pkclick import CSV, SampleList
import pklib.pkcsv as csv
import phenotool.options as OPTIONS

@click.command(no_args_is_help=True, hidden=True)
@click.argument('files', nargs=-1, type=click.File())
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--samplefiles', type=SampleList(mode='rb'), help="NOT IMPLEMENTED.") # TODO: We need a way to set the input files; especially when the samples are split across multiple files.
def main(files, columns, samplefiles):
	"""NOT IMPLEMENTED YET. Output CSV file suitable for PEP.

The Portable Encapsulated Projects (PEP for short) community effort to facilitate the portability, reusability and
durability of sample metadata.

\b
For more on the PEP community effort, please refer to:
http://pep.databio.org/en/latest/
"""
	import pkpheno as Pheno
	click.echo(samplefiles)
	assert True, "."
	pheno = Pheno.PEP(csv.DictReader(files[0]), phenovars=columns)
	for fileobj in files[1:]:
		pheno_new = Pheno.PEP(csv.DictReader(fileobj), phenovars=columns)
		pheno = pheno.combine_first(pheno_new)
	pheno.write()

#@main.command()
#def test1():
#	""""""
#	pass

#@main.command()
#def test2():
#	""""""
#	pass


