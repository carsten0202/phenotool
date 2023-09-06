
###########################################################
#
# ---%%%  Phenotool: Command 'textfile' main file  %%%---
#

###########################################################
#
# --%%  Outputting customizable text files with Phenotypes  %%--

import click
import logging
import pandas as pd
import sys

from phenotool import EPILOG, OPTIONS, Phenotype
from pklib.pkclick import CSV, isalFile, SampleList
import pklib.pkcsv as csv

logger = logging.getLogger(__name__)



#
# -%  Plain Text File Command; Main Version  %-

@click.command(no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=isalFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--csv', 'formatflag', flag_value='csv', help=OPTIONS.csv)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--tsv', 'formatflag', flag_value='tsv', default=True, help=OPTIONS.tsv)
def textfile(obj, files, columns, formatflag, samples):
    """Output phenotypes in customizable text format."""
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = TextFile(fobj, dialect=dialect, phenovars=columns, samples=samples)
        try: obj['pheno'] = obj['pheno'].combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    obj['pheno'].write(sep=formatflag)



#
# -%  Plain Text File Command; Chained Version  %-

@click.command(name="textfile", no_args_is_help=True, epilog=EPILOG.chained)
@click.pass_obj
@click.argument('files', nargs=-1, type=isalFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--csv', 'formatflag', flag_value='csv', help=OPTIONS.csv)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--tsv', 'formatflag', flag_value='tsv', default=True, help=OPTIONS.tsv)
def textfile_chain(obj, files, columns, formatflag, samples):
    """Output phenotypes in customizable text format."""
    def processor(pheno):
        if obj.get('to_be_deleted'):
            pheno.df = pheno.drop(obj['to_be_deleted'], axis='columns')
        pheno = pheno.to_textfile()
        pheno.write(sep=formatflag)
        return pheno

    try: obj['args']
    except KeyError: obj['args'] = dict()
    if columns:
        obj['args']['phenovars'] = list(dict.fromkeys(obj['args'].get('phenovars', []) + obj.get('to_be_deleted', []) + columns)) # Clever little trick to get unique list
    if samples:
        obj['args']['samples'] = list(dict.fromkeys(obj.get('samples', []) + samples))
    obj['constructor'] = obj.get('constructor', TextFile)
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = obj['constructor'](fobj, dialect=dialect, **obj['args'])
        try: obj['pheno'] = obj['pheno'].combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    return processor




###########################################################
#
# --%%  DEFINE: TestFile class  %%--

class TextFile(Phenotype):
    """Holds phenotypes in a custom text format for flexible output.
    """
    __name__ = "TextFile"
    MAGIC_COLS = {'ID'    : Phenotype.MAGIC_COLS['IID'],
                  'SEX'   : Phenotype.MAGIC_COLS['SEX']}
    mkey_id  = "ID" # Also the index, so must be unique.
    FORMATFLAGS = {'csv' : ',',
                   'tsv' : '\t'}

    def __init__(self, *args, **kwargs):
        """Init the TextFile object."""
        super().__init__(*args, **kwargs)

    def to_textfile(self):
        """No conversion necessary; self is already TextFile."""
        return self

    def write(self, sep="tsv", dest=sys.stdout, *args, **kwargs):
        """Output with support for textfile formatflags."""
        super().write(*args, sep=self.FORMATFLAGS[sep], **kwargs)



