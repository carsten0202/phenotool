
###########################################################
#
# ---%%%  Phenotool: Command 'rvtest' main file  %%%---
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
# -%  Snptest Command; Main Version  %-

@click.command(no_args_is_help=True)
@click.pass_obj
@click.argument('files', nargs=-1, type=gzFile(mode='rb'))
@click.option('-c', '--covariates', type=CSV(), default="", help=OPTIONS.covariates)
@click.option('-p', '--phenotypes', type=CSV(), default="", help=OPTIONS.phenotypes)
@click.option('-s', '--samples', type=SampleList(mode='rb'), help=OPTIONS.samples)
def snptest(obj, files, covariates, phenotypes, samples):
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
    for fobj in files:
        if (dialect := csv.sniff(fobj)) is None:
            fobj = csv.DictReader(fobj)
        pheno_new = Snptest(fobj, dialect=dialect, covariates=covariates, phenovars=phenotypes, samples=samples)
        try: obj['pheno'] = obj['pheno'].combine_first(pheno_new)
        except KeyError: obj['pheno'] = pheno_new
    obj['pheno'].write()



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
            pheno.df = pheno.drop(ctx.obj['to_be_deleted'], axis='columns')
        pheno = pheno.to_snptest(covariates = covariates)
        pheno.write()
        return pheno

    ctx.obj['phenovars'] = list(dict.fromkeys(ctx.obj.get('phenovars', []) + covariates + phenotypes)) # Clever little trick to get unique list
    if samples:
        ctx.obj['samples'] = list(dict.fromkeys(ctx.obj.get('samples', []) + samples))
    ctx.obj['constructor'] = ctx.obj.get('constructor', Snptest)
    ctx.obj['pheno'] = ctx.obj['constructor'](csv.DictReader(files[0]), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
    for fileobj in files[1:]:
        pheno_new = ctx.obj['constructor'](csv.DictReader(fileobj), phenovars=ctx.obj['phenovars'], samples=ctx.obj.get('samples'))
        ctx.obj['pheno'] = ctx.obj['pheno'].combine_first(pheno_new)
    return processor




##################################################
#
# --%%  RUN: Define CLASS Snptest  %%--

# TODO: We currently have no way to properly calculate the missing data proportion of each individual. (The 'missing' Column)
# TODO: Ok, we need to specify in detail the phenotypes and covariates. 
#   Put them in __init__? self.covar = stuff ?
#   Also means we must read files directly as Snptest....

#@pd.api.extensions.register_dataframe_accessor("pheno")
class Snptest(Phenotype):
    """Holds phenotypes in a format appropriate for Snptest
    Documentation can be found here:
    
    https://www.well.ox.ac.uk/~gav/snptest/#input_file_formats
    https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html
    """
    __name__ = "Snptest"
    MAGIC_COLS = {'ID_1'    : Phenotype.MAGIC_COLS['IID'],
                  'ID_2'    : Phenotype.MAGIC_COLS['FID'],
                  'missing' : Phenotype.MAGIC_COLS['MISSING'],
                  'sex'     : Phenotype.MAGIC_COLS['SEX']}

    mkey_id    = "ID_1" # Also the index, so must be unique.
    mkey_altid = "ID_2"
    mkey_sex   = "sex"

    def __init__(self, *args, phenovars=[], covariates=[], **kwargs):
        """Init the Snptest object. Phenovars are labelled as phenotypes unless given in covariates."""
        super().__init__(*args, phenovars = phenovars + covariates, **kwargs)
        self._obj[self.mkey_altid] = self._obj.index
        self._obj['missing'] = self._obj.get('missing', pd.Series(pd.NA * self._obj.index.size, index=self._obj.index))
        self.sex = self._obj.get(self.mkey_sex, pd.NA * self._obj.index.size)
        self.covariates = covariates
        self.phenotypes = [p for p in self.colnames_normal if p not in covariates]
        Snptest._validate(self)

    @staticmethod
    def _validate(self, warn=False):
        # We call the parent function without the meta-magic cols, for these should be validated here.
        Phenotype._validate(self[[self.mkey_sex] + self.colnames_normal], warn)

    @property
    def coltype(self):
        """Set the column type for the pseudo-type-header in sample files.
        Values are:
        0 (for the identifier column(s) and for the 'missing' column)
        D (for a column containing discrete values, e.g. a set of strings)
        P or C - for columns containing continuous value - each value must be numerical, or a missing value.
        B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.
        Sex of samples should be in column 'sex' of type D. Values in the column must be "f", "m", "male" or "female".
        """
# WARNING: No proper support for binary traits yet... must only be '0 = control', or '1 = case'
# TODO: Use the types from the arguments and input files.
        stype = pd.Series(index=self._obj.columns, dtype='object')
        for colname in self._obj.columns:
            Dtype = self._obj[colname].dropna().convert_dtypes().dtype
            if colname in self.colnames_magic:
                stype[colname] = 'D' if colname == self.mkey_sex else '0'
            elif Dtype.name == 'category':
                stype[colname] = 'B' if all(self._obj[colname].dropna().isin([0,1,'0','1','Case','case','Control','control'])) else 'D'
            elif pd.api.types.is_numeric_dtype(Dtype):
                stype[colname] = 'C'
            elif Dtype.name == 'string':
                stype[colname] = 'D'
        stype[self.covariates] = 'C'
        stype[self.phenotypes] = 'P'
        for nacol in stype.index[stype.isna()]:
            logger.warning(f"Unable to properly identify datatype for column '{nacol}'. Setting to 'D' (Discrete values).")
            stype[nacol] = 'D'
        return stype

    @property
    def covariates(self):
        return self._covariates

    @covariates.setter
    def covariates(self, value):
        value = pd.Index(value)
        self._covariates = value.intersection(self._obj.select_dtypes(include='number').columns)

    @property
    def phenotypes(self):
        return self._phenotypes

    @phenotypes.setter
    def phenotypes(self, value):
        value = pd.Index(value)
        self._phenotypes = value.intersection(self._obj.select_dtypes(include='number').columns)

    def combine_first(self, other, *args, **kwargs):
        """Super(), then set covariates/phenotypes."""
        obj = super().combine_first(other, *args, **kwargs)
        obj.covariates = self.covariates.union(other.covariates)
        obj.phenotypes = self.phenotypes.union(other.phenotypes)
        return obj

    def to_snptest(self, covariates=None):
        """Checks if covariates is not None and uses self.covariates instead."""
        if covariates:
            return super().to_snptest(covariates)
        else:
            return super().to_snptest(self.covariates)

    def write(self, *args, dest=sys.stdout, **kwargs):
# NOTE: All phenotypes should appear after the covariates in this file.
        self = self.set_columns()
        print(' '.join([self.mkey_id] + self._obj.columns.to_list()))
        print(' '.join(['0'] + self.coltype.to_list()))
        super().write(dest, *args, sep=' ', na_rep='NA', header=False, **kwargs)

# --%%  END: Define CLASS Snptest  %%--
#
##################################################




