
###########################################################
#
# ---%%%  Phenotool: Command 'textfile' main file  %%%---
#

import click

import phenotool.epilog as EPILOG
import phenotool.options as OPTIONS
from pklib.pkclick import CSV, gzFile, SampleList
import pklib.pkcsv as csv


#
# -%  Plain Text File Command; Main Version  %-

@click.group(invoke_without_command=True, no_args_is_help=True)
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--csv', 'formatflag', flag_value='csv', help=OPTIONS.csv)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--tsv', 'formatflag', flag_value='tsv', default=True, help=OPTIONS.tsv)
def textfile(files, columns, formatflag, samples):
    """Output phenotypes in customizable text format."""
    import pkpheno as Pheno
    import pandas as pd
    pheno = Pheno.TextFile(csv.DictReader(files[0]), phenovars=columns, samples=samples)
    for fileobj in files[1:]:
        pheno_new = Pheno.TextFile(csv.DictReader(fileobj), phenovars=columns, samples=samples)
        pheno = pheno.combine_first(pheno_new)
    pheno.write(sep=formatflag)



#
# -%  Plain Text File Command; Chained Version  %-

@click.command(name="textfile", no_args_is_help=True, epilog=EPILOG.chained)
@click.pass_context
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('--csv', 'formatflag', flag_value='csv', help=OPTIONS.csv)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
@click.option('--tsv', 'formatflag', flag_value='tsv', default=True, help=OPTIONS.tsv)
def textfile_chain(ctx, files, columns, formatflag, samples):
    """Output phenotypes in customizable text format."""
    def processor(pheno):
        if ctx.obj.get('to_be_deleted'):
            pheno = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
        pheno = pheno.to_textfile()
        pheno.write(sep=formatflag)
        return pheno

    ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + columns)) # Clever little trick to get unique list
    if samples:
        ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
    from pkpheno import TextFile
    ctx.obj['constructor'] = ctx.obj.get('constructor', TextFile)
    ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
    for fileobj in files[1:]:
        pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
        ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
    return processor


