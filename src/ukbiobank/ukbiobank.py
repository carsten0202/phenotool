

import itertools
import logging
import pandas as pd
import re
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

#import pklib.pkcsv as csv
from phenotool import Phenotype



##################################################
#
# --%%  CLASS: UKBioBank    %%--

class UKBioBank(Phenotype):
    """Holds phenotypes in a format appropriate for UKBioBank."""
    __name__ = "UKBioBank"
    MAGIC_COLS = {
        "EID":       ["f.eid","eid"], # In Actuality, the ID 'column' isn't a column, it's the index.
        "PAT":       Phenotype.MAGIC_COLS["PAT"],
        "MAT":       Phenotype.MAGIC_COLS["MAT"],
        "SEX":       ["f.31.0.0","31-0.0"], # UKB encodes as 0=female, but SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
    }

    mkey_id      = "EID" # Also the index, so must be unique.
    mkey_sex     = "SEX"

    def __init__(self, iterable, *args, phenovars, samples=[], sexcol=None, values=None, **kwargs):
        """
        iterable: An iterable with data...
        phenovars: The UKBiobank datafield(s) to extract. (Named for compatibility with ancestor classes).
        """
# TODO: Fix 'sex'
# NOTE: The second digit in datafields is called an 'instance'.
# NOTE: The third is the 'array index'.
        sexcol = 31 # TODO: Fix the sex col
        super().__init__(iterable, *args, usecols=lambda x: any([re.match(f'(f\D|){y}\D', x) for y in phenovars + ['ei', sexcol]]), phenovars=phenovars, samples=samples, **kwargs)

    def _conform_columns(self, columns=[]):
        """Overloads generic to set standardized names of ukb columns regardless of tab/csv origin.
        """
        self = super()._conform_columns(columns=columns)
        patone = re.compile("f\.")
        pattwo = re.compile("[-.]")
        self.df = self.df.rename(columns=lambda label: pattwo.sub("_",patone.sub('f',label)))
        logger.debug(f"UKBioBank: Renamed columns {list(self._obj.columns)}")
        return self

    @property
    def sex(self):
        """Returns the SEX in a systematic way (male/female) for querying."""
        out = pd.Series([pd.NA] * self._obj.index.size, name=self.mkey_sex, index=self._obj.index)
        out[self[self.mkey_sex].pkisin(['1'])] = "male"
        out[self[self.mkey_sex].pkisin(['0'])] = "female"
        return out

    def dc13toDate(self, fields):
        """Convert pseudo-dates in data coding 13 format to pythonic dates for specified fields.
        Return: Copy of self where 'fields' are converted."""
        def asPeriod(value):
            try:
                value = pd.NA if int(value) in [-1, -3] else pd.Period(freq='M', year=int(value), month=value %1 * 12)
            except TypeError:
                pass
            return value

        mycols = self.field2cols('20008')
        self[mycols] = self[mycols]._obj.applymap(asPeriod)
        logger.debug(f"dc13toDate: Converted subset (firstrow) = {self[mycols]._obj.iloc[0].to_list()}")
        return self

    def drop(self, labels=None, index=None, columns=None, *args, **kwargs):
        if labels is not None:
            labels = self.field2cols(labels)
            return super().drop(labels, *args, **kwargs)
        if index is not None:
            index = self.field2cols(index)
            return super().drop(index, *args, **kwargs)
        if columns is not None:
            columns = self.field2cols(columns)
            return super().drop(columns, *args, **kwargs)
        raise TypeError("Drop should have one and only one of labels, index or columns speficied.")

    def field2cols(self, fields):
        """Find all column names containing 'field' using regex."""
        if isinstance(fields, str):
            fields = [fields]
        fields = [f"^\D*{field}\D" for field in fields]
        return super().field2cols(fields)

    def findinfield(self, fields, values, instances=None, *args, mask=None):
        """Find value in indicated fields.

        fields: The field(s) to search through. Forwarded to field2cols.
        values: The value(s) to search for. Forwarded to pkisin.
        instances: A pd.Series with same index as self holding the instances to use.
        mask: An optional mask to apply before search. Forwarded to pd.mask()

        Return: A pd.Series with True/False for each row describing if row contained the value.

        Note: It seems -1 and -3 are consistently used to indicate 'missing' in UKB. This is implemented here.
        """
        # Ok, here's a serious bug. Value 'NA' in UKB doesn't mean 'pd.NA', rather it means False.
        mymask = self[self.field2cols(fields)].pkisin(['-1','-3'])
        if instances is not None:
            df = self[self.field2cols(fields)]
            mymask = mymask | ~pd.DataFrame([pd.Series(df.columns, index=df.columns).str.contains(f"_{'_|_'.join(i)}_") for i in instances], index=instances.index)
        if mask is not None:
            mymask = mymask | mask
        out = super().findinfield(fields=fields, values=values, *args, mask=mymask)
        return out

    def findinterpolated(self, field, other, values):
        """Perform lookup in UKB interpolated data

        field: The field from which data should be returned.
        other: The field where the lookup should be done.
        values: The values to find in other.

        Return: A pandas obj with values from field corresponding to values were found in other.
        """
        cols1 = self.field2cols(field)
        cols2 = self.field2cols(other)
        mask = self[cols2].pkisin(values)
        mask = mask.rename(columns=dict(zip(cols2,cols1)))
        out = self._obj[cols1].mask(~mask, other=pd.NA)
        out = out.mask(out.isin([-1,-3,'-1','-3']), other=pd.NA).dropna(axis='columns', how='all')
        logger.debug(f"findinterpolated: field={field}; other={other}; values={values}; return={out.to_dict()}")
        return out

#    def getfield(self, fields):
#        """Lookup """
#        self

