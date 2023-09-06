
###########################################################
#
# ---%%%  Phenotype: Handling Phenotype information  %%%---
#

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

import copy
import logging
from numbers import Number
import numpy as np
import pandas as pd
import sys
import warnings

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import pklib

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  RUN: Define CLASS pheno  %%--

class Phenotype:
    """A handy-dandy object for storing phenotypes. Will be cool to also import into other scripts when working with phenotypes"""
    __name__ = "Phenotype"
    MAGIC_COLS = {"IID":     ["ID","id","ID_1","Id_1","id_1","IID","Iid","iid","ParticID","Particid","particid","sample"],
                  "FID":     ["FID","Fid","fid","ID_2","Id_2","id_2"],
                  "MAT":     ["MAT"],
                  "MISSING": ["MISSING","Missing","missing"],
                  "PAT":     ["PAT"],
                  "SEX":     ["SEX","Sex","sex","GENDER","Gender","gender"], # SEX should be encoded with males as '1' or 'M' or 'male'; females as '2' or 'F' or 'female'
                  "SID":     ["SID"]}

    mkey_id    = "IID" # Also the index, so must be unique.
    mkey_sex   = "SEX"

    def __init__(self, iterable, *args, phenovars=[], samples=[], **kwargs):
        """
        iterable:  Someting iterable. Possibly a 'Phenotype' class object.
        phenovars: A list of variables to output. If empty it must default to all columns in 'iterable'.
        """
        try:
            self.df = pd.read_csv(iterable, *args, **kwargs)
        except:
            try: self.df = pd.DataFrame(iterable)
            except Exception as ex:
                logger.critical(f'Phenotype: Unable to parse input for Phenotype constructor. Iterable expected, received {iterable}')
                sys.exit(ex)
        logger.info(f"{self.__name__}: Parsing file with columns[:10] = {self.colnames[:10]}")
        logger.debug(f"{self.__name__}: Parsing file with columns: {self.colnames}")
        logger.info(f"{self.__name__}: Searching for vars = {phenovars}")
        if samples:
            self.samples = samples
        self = self._set_magic_kcol()
        self.df = self.df.set_index(self.mkey_id)
        self = self._conform_columns(columns=self.field2cols(phenovars))
        Phenotype._validate(self)
        logger.debug(f"{self.__name__}: Columns after __init__ = {self.columns.to_list()}")

# This sample data shit should probably be moved to a class which deals with sample data. It shouldn't be here in the generic...
        if self.is_sample_data():
            self._sampletypes = self.index[0]
            self.df = self.drop(self.index[0])

    def _conform_columns(self, columns=[]):
        """Return: All magic columns in input + specified columns in that order. Also converts data types and NAs.
        columns: a list of column names to extract. Default: All.
        """
        columns = self.colnames_translate(columns) if columns else self.colnames
        index = pd.Index(self.colnames_magic).append(pd.Index(columns)).drop_duplicates().to_list()
        cols = [col for col in index if col in self.df]
        self.df = self.df.loc[:,cols]
        self.df.replace('NA', np.NaN, inplace=True)
        for col in self.df.select_dtypes('object').columns:
            self._obj[col] = pd.to_numeric(self._obj[col], errors='ignore')
        self.df = self.df.convert_dtypes()
        return self

    def __getitem__(self, key):
        """Redirects index operations to the _obj attribute."""
        out = copy.deepcopy(self)
        out._obj = out._obj[key]
        return out

    def __repr__(self):
        """Use the _obj's __repr__ function instead; ie show the pandas object."""
        return self._obj.__repr__()


    def __setitem__(self, key, value):
        """Redirects index-based assignment operations to the _obj attribute."""
        self._obj[key] = value

    def _set_magic_kcol(self):
        for k, v in self.MAGIC_COLS.items():
            self.df = self.df.rename(dict(zip(v, [k] * len(v))), axis=1)
        if self.mkey_id not in self._obj:
            self.df = self.df.rename(columns={ self.columns[0]: self.mkey_id })
            logger.warning(f"Using first column as primary identifier as no column was found with recognized headers: {', '.join(self.MAGIC_COLS[self.mkey_id])}.")
        assert self.mkey_id in self.df, f"Unable to identify primary identifier in one or more input files. Please make sure that each input file has a column with one of these recognized headings: {', '.join(self.MAGIC_COLS[self.mkey_id])}"
        return self

    @staticmethod
    def _validate(self, warn=False):
        # Check What? 
        # We should check the SEX column, if it is there... And we need to estimate 0,1,2 errors (when 1 cn mean both males (psam) and females (rvtest)).
        # We should check that all columns have headers and that the headers are strings (numerical headers usually indicates a problem).
        # We should check that there is a header (difficult...)
        dup = self.df.index.duplicated()
        assert not dup.any(), f"One or more input files had duplicated primary identifiers ({self._obj.index[dup].unique().to_list()}). The primary identifiers must be unique within each input file."
        if warn:
            samples_notok = pd.Series(pd.NA * self.index.size, index=self.index)
            columns_notok = pd.Series(pd.NA * len(self.columns), index=self.columns)
            for name,row in self.df.iterrows():
                if row.isna().all():
                    samples_notok[name] = name
            samples_notok = samples_notok.dropna().to_list()
            for name,col in self.df.items():
                columns_notok[name] = name if col.isna().all() else pd.NA
            for i,sample in enumerate(samples_notok):
                if i == 0:
                    first = samples_notok[0:3] + ['...'] if len(samples_notok) > 3 else samples_notok
                    logger.warning(f"Dataset contains {len(samples_notok)} sample(s) without any associated data {first}.")
                logger.info(f"Sample '{sample}' has no associated data; all it's variables set to missing values.")
            columns_notok = columns_notok.dropna().to_list()
            for i,column in enumerate(columns_notok):
                if i == 0:
                    first = columns_notok[0:3] + ['...'] if len(columns_notok) > 3 else columns_notok
                    logger.warning(f"Dataset contains {len(columns_notok)} phenotype(s) without any associated data {first}.")
                logger.info(f"Phenotype '{column}' has no data associated with any subjects present in the data.")

    @property
    def colnames(self):
        """Return: All colnames in _obj as list."""
        return self._obj.columns.to_list()

    @property
    def colnames_magic(self):
        """Return: All magic colnames supported by class as list."""
        return list(self.MAGIC_COLS.keys())

    @property
    def colnames_normal(self):
        """Return: All colnames in _obj which are not magical as list."""
        return [col for col in self.columns if col not in self.colnames_magic]

    @property
    def columns(self):
        """Returns the columns of _obj."""
        return self._obj.columns

    @property
    def df(self):
        """Return the pheno DataFrame."""
        return self._obj

    @df.setter
    def df(self, value):
        """Setter: Set the df (_obj) to value."""
        try: self._obj = value
        except Exception as ex:
            sys.exit(f"{ex} in df.setter...")

    @property
    def index(self):
        """Return the index of _obj."""
        return self._obj.index

    @property
    def samples(self):
        """Synonym for self.index."""
        return self.index

    @samples.setter
    def samples(self, value):
        """Return: All specified samples in specified order."""
        try: self.df = self.df.reindex(value)
        except Exception as ex:
            sys.exit(f"{ex} in samples.setter...")

    @property
    def sex(self):
        """Returns the SEX in a systematic way (male/female or child-dependent) for querying."""
        out = pd.Series(np.zeros(self._obj.index.size) + np.nan, name=self.mkey_sex, index=self._obj.index)
        out[self[self.mkey_sex].pkisin(['1', 'M', "male"])] = "male"
        out[self[self.mkey_sex].pkisin(['2', 'F', "female"])] = "female"
        return out

    @sex.setter
    def sex(self, value):
        """Sets the sex column in a somewhat error-robust way."""
        val = pd.Series(value, index=value.index if isinstance(value, pd.Series) else self.index, dtype="category")
        val = val.cat.rename_categories(dict(zip([2, 'F', 'f', 'FEMALE', 1, 'M', 'm', 'MALE'], ['female']*4 + ['male']*4)))
        self._obj[self.mkey_sex] = val

    def colnames_translate(self, columns):
        """Translate names in columns to the canonical magic column names used by self."""
        outcols = []
        for col in columns:
            mcol = [k for k in self.colnames_magic if col in self.MAGIC_COLS[k]]
            outcols.extend(mcol if mcol else [col])
        return outcols

    def combine_first(self, other):
        """Return: Combined DataFrame from self and other."""
        if isinstance(other, Phenotype):
            self.df = self.df.combine_first(other=other.df)
        else:
            self.df = self.df.combine_first(other=other)
        self.df = self.df.fillna(np.NaN)
        return self

    def derive_minmax(self, columns=None):
        """Return: DataFrame with results scaling according to algorithm and cols given by 'columns'.
        Default: All normal columns."""
        columns = self.colnames_translate(columns) if columns else self.colnames_normal
        df = pd.DataFrame()
        for col in columns:
            df[col] = pklib.scaler_min_max(self.df[col])
        return df

    def derive_rankINT(self, columns=None):
        """Return: DataFrame with results from a rank-based inverse normally transformation on cols given by 'columns'.
        Default: All normal columns."""
        columns = self.colnames_translate(columns) if columns else self.colnames_normal
        df = pd.DataFrame()
        for col in columns:
            df[col] = pklib.rank_int(self.df[col])
        return df

    def derive_zscores(self, columns=None):
        """Return: DataFrame with results scaling according to algorithm and cols given by 'columns'.
        Default: All normal columns."""
        columns = self.colnames_translate(columns) if columns else self.colnames_normal
        df = pd.DataFrame()
        for col in columns:
            df[col] = pklib.scaler_min_max(self.df[col])
        return df

    def drop(self, *args, **kwargs):
        """NOTE: I changed this guy from dropping inplace. Inplace was stupid as you could simply add argument 'inplace=True' to get it."""
        return self.df.drop(*args, **kwargs)

    def field2cols(self, fields):
        """Find all column names containing 'field' using regex."""
        out = pd.Series()
        cols = pd.Series(self.columns)
        if not isinstance(fields, list):
            fields = [fields]
        for field in fields:
            out = pd.concat([out, cols[cols.str.contains(str(field), regex=True)]])
        logger.debug(f"field2cols: Fields={fields}; Found cols={out.to_list()}.")
        return out.to_list()

    def findfield(self, fields):
        """Find all column names containing fields and return those as DataFrame."""
        out = self.df[self.field2cols(fields)]
        return out

    def findinfield(self, fields, values, *args, mask=None):
        """Find values in indicated fields.

        fields: The field(s) to search through. Forwarded to field2cols.
        values: The value(s) to search for. Forwarded to pkisin.
        mask: An optional mask to apply before search. Forwarded to pd.mask()

        Return: A pd.Series with True/False for each row describing if row contained the value.
        """
        out = pd.Series(pd.NA, index=self.index, dtype='boolean')
        cols = self.field2cols(fields)
        if self.columns.isin(cols).any():
            df = self[cols]
            if mask is None:
                mask = pd.DataFrame(False, index=df.index, columns=df.columns)
            else:
                logger.debug(f"findinfield: Mask = {mask.to_dict()}")
            out[df.index] = df.pkisin(values).any(axis='columns')
            out[df._obj.mask(mask, pd.NA).isna().all(axis='columns')] = pd.NA
        logger.debug(f"findinfield: Scanning for values={values}; Found = {out.to_dict()}.")
        return out

    def is_sample_data(self):
        """0 (for the first identifier column)
        D (for a column containing discrete values, e.g. a set of strings)
        P or C - for columns containing continuous value - each value must be numerical, or a missing value.
        B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.
        """
        list_ = [0,"0","B","C","D","P"]
        return all(self._obj.iloc[0].isin(list_))

    def pkisin(self, columns, values):
        """Convienience function to check for 'values' using both numeric and string in df.isin()."""
        if isinstance(values, str):
            values = [values]
        out = values
        for value in values:
            if isinstance(value, str) and value.isnumeric():
                out.append(float(value))
            elif isinstance(value, Number):
                out.append(str(value))
        return self.df[columns].isin(out)

    def to_psam(self):
        """Convert Phenotype Class to Class Psam for Plink output."""
        from .plink import Psam
        obj = self.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None), getattr(self, 'mkey_mat', None), getattr(self, 'mkey_pat', None)], errors='ignore')
        obj[Psam.mkey_id] = self.index
        obj[Psam.mkey_sex] = self.sex
        try: obj[Psam.mkey_altid] = self._obj[self.mkey_altid]
        except AttributeError:
            pass
        try: obj[Psam.mkey_mat] = self._obj[self.mkey_mat]
        except AttributeError:
            pass
        try: obj[Psam.mkey_pat] = self._obj[self.mkey_pat]
        except AttributeError:
            pass
        psam = Psam(obj)
        return psam

    def to_rvtest(self):
        """Convert Phenotype Class to Class RVtest for RVtest output."""
        from .rvtest import RVtest
        obj = self._obj.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None), getattr(self, 'mkey_mat', None), getattr(self, 'mkey_pat', None)], errors='ignore')
        obj[RVtest.mkey_id] = self.index
        obj[RVtest.mkey_sex] = self.sex
        try: obj[RVtest.mkey_altid] = self._obj[self.mkey_altid]
        except AttributeError:
            pass
        try: obj[RVtest.mkey_mat] = self._obj[self.mkey_mat]
        except AttributeError:
            pass
        try: obj[RVtest.mkey_pat] = self._obj[self.mkey_pat]
        except AttributeError:
            pass
        rvtest = RVtest(obj)
        return rvtest

    def to_snptest(self, covariates=[]):
        """Convert Phenotype Class to Class Snptest for sample file output."""
        from .snptest import Snptest
        obj = self.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None)], errors='ignore')
        obj[Snptest.mkey_id] = self.index
        obj[Snptest.mkey_sex] = self.sex
        try: obj[Snptest.mkey_altid] = self._obj[self.mkey_altid]
        except AttributeError:
            pass
        try: snptest = Snptest(obj, covariates=covariates)
        except TypeError:
            logger.error("to_snptest: Something went wrong in converting to Snptest. Did you set your covariates correctly?")
            raise
        return snptest

    def to_textfile(self):
        """Convert Phenotype Class to Class TextFile for custom file output."""
        from .textfile import TextFile
        obj = self.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None)], errors='ignore')
        obj[TextFile.mkey_id] = self.index
        obj[TextFile.mkey_sex] = self.sex
        try: obj[TextFile.mkey_altid] = self._obj[self.mkey_altid]
        except AttributeError:
            pass
        textfile = TextFile(obj)
        return textfile

    def write(self, dest=sys.stdout, *args, **kwargs):
        """Output self._obj to dest."""
        self._validate(self, warn=True)
        self.df.to_csv(dest, *args, **kwargs)

# --%%  END: Define CLASS pheno  %%--
#
##################################################




##################################################
#
# --%%  RUN: Define CLASS PEP for Portable Encapsulated Projects  %%--

class PEP(Phenotype):
    """Portable Encapsulated Projects

    """
    __name__ = "PEP"
    MAGIC_COLS = {
        "sample_name":     Phenotype.MAGIC_COLS["IID"] + ["sample_name"], # In Actuality, the ID 'column' isn't a column, it's the index.
    }
    mkey_id    = "sample_name" # Also the index, so must be unique.

    def __init__(self, *args, **kwargs):
        """"""
        super().__init__(*args, **kwargs)
        PEP._validate(self)

    def write(self, *args, quoting=pklib.QUOTE_ALL, **kwargs):
        super().write(*args, quoting=quoting, **kwargs)

# --%% END: Defint CLASS PEP for Portable Encapsulated Projects  %%--
#
##################################################


##################################################
#
# --%%  RUN: Other Functions & Constructors  %%--

def rank_INT(series, c=3.0/8, stochastic=True):
    """Perform rank-based inverse normal transformation on pandas series.

    If stochastic is True ties are given rank randomly, otherwise ties will
    share the same value. NaN values are ignored.
    Args:
        param1 (pandas.Series):   Series of values to transform
        param2 (Optional[float]): Constand parameter (Bloms constant)
        param3 (Optional[bool]):  Whether to randomise rank of ties
    Returns:
        pandas.Series
    Inspired by:
        https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html
    """
    import pandas as pd
    import scipy.stats as ss
    def rank_to_normal(rank, c, n):
        x = (rank - c) / (n - 2*c + 1) # Standard quantile function
        return ss.norm.ppf(x)

    # Check input
    assert isinstance(series, pd.Series)
    assert isinstance(c, float)

    np.random.seed(123) # Set seed
    orig_idx = series.index # Take original series indexes
    series = series.loc[~pd.isnull(series)] # Drop NaNs
    if stochastic:
        series = series.loc[np.random.permutation(series.index)] # Shuffle by index
        rank = ss.rankdata(series, method="ordinal") # Get rank, ties are determined by their position in the series (hence why we randomised the series)
    else:
        rank = ss.rankdata(series, method="average") # Get rank, ties are averaged
    rank = pd.Series(rank, index=series.index) # Convert numpy array back to series
    transformed = rank.apply(rank_to_normal, c=c, n=len(rank)) # Convert rank to normal distribution
    return transformed[orig_idx]

# --%% END: Other Functions & Constructors  %%--
#
##################################################

