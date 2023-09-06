

import itertools
import logging
import pandas as pd
import re
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

from phenotool import Phenotype



##################################################
#
# --%%  CLASS: UKBioBank    %%--

class UKBioBank(Phenotype):
    """Holds phenotypes in a format appropriate for UKBioBank."""
    __name__ = "UKBioBank"
    MAGIC_COLS = {
        "EID":       ["f.eid","eid"], # In Actuality, the ID 'column' isn't a column, it's the index.
        "SEX":       ["f.31.0.0","31-0.0"], # UKB encodes as 0=female, but SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
    }

    mkey_id      = "EID" # Also the index, so must be unique.
    mkey_sex     = "SEX"

    sex_dict = {
        "genetic":  ["f.22001.0.0","22001-0.0"],
        "registry": ["f.31.0.0","31-0.0"],
    }

    def __init__(self, iterable, *args, instances=[], phenovars, samples=[], sexcol=None, values=None, **kwargs):
        """
        iterable: An iterable with data...
        phenovars: The UKBiobank datafield(s) to extract. (Named for compatibility with ancestor classes).
        """
# NOTE: The second digit in datafields is called an 'instance'.
# NOTE: The third is the 'array index'.
        self.MAGIC_COLS[self.mkey_sex] = self.sex_dict.get(sexcol, [])
        if instances:
            phenovars = [f"{p}\D[{''.join(instances)}]" for p in phenovars]
        col_fun = lambda x: any([re.match(f'(f\D|){p}\D', x) for p in phenovars] + [x in list(itertools.chain.from_iterable(self.MAGIC_COLS.values()))])
        super().__init__(iterable, *args, usecols=col_fun, phenovars=phenovars, samples=samples, **kwargs)

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
        try:
            out[self.pkisin(self.mkey_sex, ['0', 'female'])] = "female"
            out[self.pkisin(self.mkey_sex, ['1', 'male'])] = "male"
        except KeyError:
            pass
        return out

    @sex.setter
    def sex(self, value):
        """Sets the sex column in a somewhat error-robust way."""
        try: value = value.squeeze()
        except: pass
        val = pd.Series(value, index=value.index if isinstance(value, pd.Series) else self.index, dtype="category")
        val = val.cat.rename_categories(dict(zip([0, 2, 'F', 'f', 'FEMALE', 1, 'M', 'm', 'MALE'], ['female']*5 + ['male']*4)))
        self.df[self.mkey_sex] = val

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
            return super().drop(labels=labels, *args, **kwargs)
        if index is not None:
            index = self.field2cols(index)
            return super().drop(index=index, *args, **kwargs)
        if columns is not None:
            columns = self.field2cols(columns)
            return super().drop(columns=columns, *args, **kwargs)
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

