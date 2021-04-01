#!/usr/bin/python

import click

from pklib.pkclick import CSV, SampleList
import pklib.pkcsv as csv
import phenotool.options as OPTIONS

@click.command(no_args_is_help=True, hidden=True)
@click.argument('files', nargs=-1, type=click.File())
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
def main(files, columns):
	"""NOT IMPLEMENTED YET. Output CSV file suitable for PEP.

The Portable Encapsulated Projects (PEP for short) community effort to facilitate the portability, reusability and
durability of sample metadata.

\b
For more on the PEP community effort, please refer to:
http://pep.databio.org/en/latest/
"""
	import pkpheno as Pheno
	assert True, "."
	pheno = Pheno.PEP(csv.DictReader(files[0]), columns=columns)
	for fileobj in files[1:]:
		pheno_new = Pheno.PEP(csv.DictReader(fileobj), columns=columns)
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


