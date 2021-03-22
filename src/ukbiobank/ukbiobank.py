#!/usr/bin/python

import click

@click.group()
def ukbiobank():
	"""UNDER DEVELOPMENT - DO NOT USE."""
	click.echo('Hello UKB World!')

@ukbiobank.command()
def test1():
	""""""
	pass

@ukbiobank.command()
def test2():
	""""""
	pass

