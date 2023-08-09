
###########################################################
#
# ---%%%  Phenotool: Command 'plink' main file  %%%---
#

import click
import logging
import pandas as pd
import sys

from phenotool import OPTIONS, Phenotype
from pklib.pkclick import CSV, gzFile, SampleList
import pklib.pkcsv as csv

logger = logging.getLogger(__name__)



#
# -%  Plink Command; Main Version  %-

@click.command(no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-f', '--fam', type=str, metavar='COLUMN', default=None, help=OPTIONS.fam)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def plink(obj, files, columns, fam, samples):
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
        if columns:
            logger.warning(f"Note that the '--columns' option is ignored when setting '--fam'.")
        columns = [fam]
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = Psam(fobj, dialect=dialect, phenovars=columns, samples=samples)
        try: obj['pheno'] = obj['pheno'].pheno.combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    obj['pheno'].write(header = False if fam else True)



#
# -%  Plink Command; Chained Version  %-

@click.command(name="plink", no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--columns', type=CSV(), default="", help=OPTIONS.columns)
@click.option('-f', '--fam', type=str, metavar='COLUMN', default=None, help=OPTIONS.fam)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def plink_chain(obj, files, columns, fam, samples):
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
    def processor(pheno):
        if obj.get('to_be_deleted'):
            pheno.df = pheno.drop(obj['to_be_deleted'], axis='columns')
        pheno = pheno.to_psam()
        pheno.write(header = False if obj.get('fam') else True)
        return pheno

    assert sum([1 for x in [columns,fam] if x]) <= 1, "'--columns' and '--fam' are mutually exclusive; please only specify one of them."
    if fam:
        if columns:
            logger.warning(f"Note that the '--columns' option is ignored when setting '--fam'.")
        columns = [fam]
        obj['fam'] = True
    try: obj['args']
    except KeyError: obj['args'] = dict()
    obj['args']['phenovars'] = list(dict.fromkeys(obj['args'].get('phenovars', []) + columns)) # Clever little trick to get unique list
    if samples:
        obj['samples'] = list(dict.fromkeys(obj.get('samples', []) + samples)) if obj.get('samples') else samples
    obj['constructor'] = obj.get('constructor', Psam)
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = obj['constructor'](fobj, dialect=dialect, **obj['args'])
        try: obj['pheno'] = obj['pheno'].combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    return processor




##################################################
#
# --%%  RUN: Define CLASS Psam (Plink2)  %%--

class Psam(Phenotype):
    """Holds phenotypes in a format appropriate for Plink2.
    Documentation can be found here:
    https://www.cog-genomics.org/plink/2.0/formats#psam
    """
    __name__ = "Psam"
    MAGIC_COLS = {
        "FID":     Phenotype.MAGIC_COLS["FID"],
        "IID":     Phenotype.MAGIC_COLS["IID"], # In Actuality, the ID 'column' isn't a column, it's the index.
        "SID":     Phenotype.MAGIC_COLS['SID'], # SID is not really supported yet
        "PAT":     Phenotype.MAGIC_COLS["PAT"],
        "MAT":     Phenotype.MAGIC_COLS["MAT"],
        "SEX":     Phenotype.MAGIC_COLS['SEX'], # SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
    }

    mkey_altid = "FID"
    mkey_id    = "IID" # Also the index, so must be unique.
    mkey_mat   = "MAT"
    mkey_pat   = "PAT"
    mkey_sex   = "SEX"

    def __init__(self, *args, **kwargs):
        """Create Phenotype object and then set the special plink columns."""
        super().__init__(*args, **kwargs)
        self._obj[self.mkey_id]    = self._obj.index
        self._obj[self.mkey_altid] = self._obj.index
        self._obj[self.mkey_pat]   = self._obj.get(self.mkey_pat, 0)
        self._obj[self.mkey_mat]   = self._obj.get(self.mkey_mat, 0)
        self._obj[self.mkey_sex]   = pd.Categorical(self._obj.get(self.mkey_sex, ['']*len(self._obj.index)))
        self._obj[self.mkey_sex]   = self._obj[self.mkey_sex].cat.rename_categories({0 : pd.NA, 'M' : 1, 'Male' : 1, 'male' : 1, 'F' : 2, 'Female' : 2, 'female' : 2})
        cols = list(dict.fromkeys([self.mkey_altid, self.mkey_id, self.mkey_mat, self.mkey_pat, self.mkey_sex] + self.colnames))
        self.df = self.df.loc[:,cols]
        Psam._validate(self)

    def to_psam(self):
        """No conversion necessary; self is already Psam."""
        return self

    def write(self, *args, dest=sys.stdout, header=True, **kwargs):
        """Like super() but handles *.fam files through header=False if self._obj only has one phenotype."""
        if header:
            print("#", end="")
        super().write(dest, *args, sep='\t', na_rep='NA', header=header, index=False, **kwargs)

# --%%  END: Define CLASS Psam (Plink2)  %%--
#
##################################################

