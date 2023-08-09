
###########################################################
#
# ---%%%  Phenotool: Command 'rvtest' main file  %%%---
#

import click
import logging
import sys

from phenotool import OPTIONS, Psam
from pklib.pkclick import CSV, gzFile, SampleList
import pklib.pkcsv as csv

logger = logging.getLogger(__name__)



#
# -%  RVtest Command; Main Version  %-
    
@click.command(no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def rvtest(obj, files, columns, samples):
    """UNTESTED; Output phenotypes in psam-like format for RVtest.

RVtest phenotype files are very similar to the psam format. They are essentially plink2 files with a few caveats. The
names 'fatid' and 'matid' are used for paternal and maternal ids and 'sex' is encoded 0=males,1=females.

\b
For more on RVtest phenotype files, please refer to:
http://zhanxw.github.io/rvtests/#phenotype-file
"""
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = RVtest(fobj, dialect=dialect, phenovars=columns, samples=samples)
        try: obj['pheno'] = obj['pheno'].pheno.combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    obj['pheno'].write()


#
# -%  RVtest Command; Chained Version  %-

@click.command(name="rvtest", no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def rvtest_chain(obj, files, columns, samples):
    """UNTESTED; Output phenotypes in psam-like format for RVtest.

RVtest phenotype files are very similar to the psam format. They are essentially plink2 files with a few caveats. The
names 'fatid' and 'matid' are used for paternal and maternal ids and 'sex' is encoded 0=males,1=females.

\b
For more on RVtest phenotype files, please refer to:
http://zhanxw.github.io/rvtests/#phenotype-file
"""
    def processor(pheno):
        if obj.get('to_be_deleted'):
            pheno.df = pheno.drop(obj['to_be_deleted'], axis='columns')
        pheno = pheno.to_rvtest()
        pheno.write()
        return pheno

    try: obj['args']
    except KeyError: obj['args'] = dict()
    obj['args']['phenovars'] = list(dict.fromkeys(obj['args'].get('phenovars', []) + columns)) # Clever little trick to get unique list
    if samples:
        obj['samples'] = list(dict.fromkeys(obj.get('samples', []) + samples))
    obj['constructor'] = obj.get('constructor', RVtest)
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = obj['constructor'](fobj, dialect=dialect, **obj['args'])
        try: obj['pheno'] = obj['pheno'].combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    return processor




##################################################
#
# --%%  RUN: Define CLASS RVtest  %%--

class RVtest(Psam):
    """Holds RVtest phenotype data
    Documentation can be found here:
    http://zhanxw.github.io/rvtests/#phenotype-file
    """
    __name__ = "RVtest"
    MAGIC_COLS = {
        "fid":     Psam.MAGIC_COLS["FID"],
        "iid":     Psam.MAGIC_COLS["IID"], # In Actuality, the ID 'column' isn't a column, it's the index.
        "fatid":   Psam.MAGIC_COLS["MAT"],
        "matid":   Psam.MAGIC_COLS["PAT"],
        "sex":     Psam.MAGIC_COLS['SEX'], # SEX should be encoded with males as '0', females as '1' and missing as 'NA'
    }

    mkey_altid = "fid"
    mkey_id    = "iid" # Also the index, so must be unique.
    mkey_mat   = "matid"
    mkey_pat   = "fatid"
    mkey_sex   = "sex"

    def __init__(self, *args, **kwargs):
        """Super() then re-encode 'sex' column: Males => 0, Females => 1."""
        super().__init__(*args, **kwargs)
        self._obj[self.mkey_sex] = self._obj[self.mkey_sex].cat.rename_categories({1 : 0, 2 : 1})
        RVtest._validate(self)

    def to_rvtest(self):
        """No conversion necessary; self is already RVtest."""
        return self

    def write(self, *args, dest=sys.stdout, **kwargs):
        """Super() skips over Psam to avoid the initial '#'."""
#        self = self.set_columns()
        super(Psam,RVtest).write(self, dest, *args, sep='\t', na_rep='NA', index=False, **kwargs)

# --%%  END: Define CLASS RVtest  %%--
#
##################################################


