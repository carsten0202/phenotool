
###########################################################
#
# ---%%%  Phenotool: Derive main file  %%%---
#

# TODO:
#   Add TextFile as default output.

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

import click
import logging
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

#import eastwood.cli as Eastwood
from phenotool import EPILOG, OPTIONS, plink_chain, rvtest_chain, snptest_chain, textfile_chain
from pklib.pkclick import CSV

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  RUN: Commands  %%--

@click.group(chain=True, no_args_is_help=True, epilog=EPILOG.chained)
@click.pass_obj
def derive(obj):
    """Derive new columns from values in existing ones.

This command group provides access to some of the mathematical functions in the pandas library.
"""

    try: obj['args']
    except KeyError: obj['args'] = dict()
    try: obj['args']['phenovars']
    except KeyError: obj['args']['phenovars'] = list()

@derive.result_callback()
@click.pass_obj
def derive_pipeline(obj, processors):
    logging.debug(f"Pipeline: Cols to be deleted: {obj.get('to_be_deleted')}")
    pheno = obj['pheno']
    for processor in processors:
        pheno = processor(pheno)

#    try: ctx.obj['args']['nrows'] = nrows
#    except KeyError:
#        ctx.obj['args'] = {'nrows': nrows}
#    ctx.obj['args']['samples'] = samples
#    ctx.obj['args']['sexcol'] = sex
#    ctx.obj['constructor'] = UKBioBank

#
# -%  Chained Output Commands  %-

# Plink output command (Chained version)
derive.add_command(plink_chain)

# RVtest output command (Chained version)
derive.add_command(rvtest_chain)

# Snptest output command (Chained version)
derive.add_command(snptest_chain)

# Textfile output command (Chained version)
derive.add_command(textfile_chain)

OPTIONS.columns_derive = "Comma separated list of columns to perform operation on."


#
# -%  DERIVE: RankINV Command  %-

@derive.command(name='rankinv', no_args_is_help=True)
@click.pass_obj
@click.option('-c', '--columns', type=CSV(), default=[], help=OPTIONS.columns_derive)
@click.option('-p', '--prefix', default="rankinv_", show_default=True, help=OPTIONS.columnprefix)
def rankinv_chain(obj, columns, prefix):
    """Perform Rank-Based Inverse Normal Transformations."""
    def processor(pheno):
        pheno[[f"{prefix}{col}" for col in columns]] = pheno.derive_rankINT(columns=columns)
        return pheno

    for col in filter(lambda col: col not in obj['args']['phenovars'], columns):
        obj['to_be_deleted'] = obj.get('to_be_deleted', list()) + [col]
    return processor



#
# -%  DERIVE: Scale Command  %-

@derive.command(name="scaling", no_args_is_help=True)
@click.pass_obj
@click.option('-c', '--columns', type=CSV(), default=[], help=OPTIONS.columns_derive)
@click.option('-p', '--prefix', default="scaled_", show_default=True, help=OPTIONS.columnprefix)
@click.option('--absmax', 'scaler', flag_value='derive_absmax', help='Maximum absolute scaling method rescales each feature to be a value between -1 and 1.')
@click.option('--minmax', 'scaler', default=True, flag_value='derive_minmax', help='Min-max feature scaling (normalization) rescales the dataset feature to a range of 0 - 1. This is default.')
@click.option('--zscores', 'scaler', flag_value='derive_zscores', help='Transforms the data into z-scores (standardization). Data becomes a distribution of values with mean 0 and standard deviation 1.')
def scaling_chain(obj, columns, prefix, scaler):
    """Perform scaling of numerical values.
    
    Implements a number of standardized numerical scaling algorithms which can be applied to one or more columns."""
    def processor(pheno):
        pheno[[f"{prefix}{col}" for col in columns]] = getattr(pheno, scaler)(columns=columns)
        return pheno

    for col in columns:
        if col not in obj['args']['phenovars']:
            obj['to_be_deleted'] = obj.get('to_be_deleted', list()) + [col]
    return processor


