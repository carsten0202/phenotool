###########################################################
#
# ---%%%  ukbiobank.py: Home of the UKBioBank Class  %%%---
#

##################################################
#
# --%%  RUN: Perform Basic Setup  %%--

import copy
import logging
import pandas as pd
import re
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import pklib.pkcsv as csv
from pkpheno import Phenotype, Psam

# --%%  END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%  CLASS: UKBioBank  %%--

class UKBioBank(Phenotype):
	"""Holds phenotypes in a format appropriate for UKBioBank."""
	__name__ = "UKBioBank"
	MAGIC_COLS = {
		"EID":     ["f.eid","eid"], # In Actuality, the ID 'column' isn't a column, it's the index.
		"PAT":     Phenotype.MAGIC_COLS["PAT"],
		"MAT":     Phenotype.MAGIC_COLS["MAT"],
		"31" :     ["31_0_0"], # This needs to be the datafield number, not the actual column name in the input
		                       # UKB encodes as 0=female, but SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
	}

	mkey_id    = "EID" # Also the index, so must be unique.
	mkey_sex   = "31"

	def __init__(self, *args, phenovars=[], **kwargs):
		""" Create an instance of class UKBioBank
		phenovars: The UKBiobank datafield(s) to extract. (Named for compatibility with ancestor classes).
		"""
# NOTE: The second digit in datafields is called an 'instance'.
# NOTE: The third is the 'array index'.
		super().__init__(*args, **kwargs)
		indexframe = pd.DataFrame([coln.split("_") for coln in self._obj.columns])
		self._obj.columns = pd.MultiIndex.from_frame(indexframe, names=['datafield','instance','array_index'])
		if fields := set(phenovars) - set(self.datafields):
			logger.warning(f"UKBioBank: One or more datafields did not exist in input data: {fields}")
		logger.debug(f"UKBioBank: Processed dataframe:\n{self._obj}")

	@classmethod
	def _read_csv_callback(cls, *args, names=None, usecols=None, **kwargs):
		"""Callback function which can be overloaded to customize the read_csv call further."""
		assert usecols, "You must use '--datafields' option to specify at least one datafield for extraction."
		magicols = '|'.join(cls.MAGIC_COLS[cls.mkey_id] + cls.MAGIC_COLS[cls.mkey_sex])
		pattern = re.compile(f"({'_|'.join(usecols)}_|{magicols})")
		usecols = lambda x: pattern.match(x)
		names = [re.sub("[f.]*([0-9]+).([0-9]+).([0-9]+)","\\1_\\2_\\3", n) for n in names]
		return pd.read_csv(*args, names=names, usecols=usecols, **kwargs)

	@property
	def datafields(self):
		return self._obj.columns.get_level_values('datafield').to_list()

	@property
	def sex(self):
		"""Returns the SEX in a systematic way (male/female) for querying.
		Overloaded for MultiIndexed columns."""
		out = pd.Series([pd.NA] * self._obj.index.size, name=self.mkey_sex, index=self._obj.index)
		out[self[self.mkey_sex,'0','0'].pkisin(['1'])] = "male"
		out[self[self.mkey_sex,'0','0'].pkisin(['0'])] = "female"
		return out

	def __getitem__(self, key):
		"""Overloaded for MultiIndexed columns."""
		out = copy.deepcopy(self)
		out._obj = out._obj.loc[:,tuple(key)]
		return out

	def _set_magic_kcol(self, *args, **kwargs):
		"""Local function to handle specific UKB problems like the multiindex."""
		df = super()._set_magic_kcol(*args,**kwargs)
		df.rename({self.mkey_sex:f'{self.mkey_sex}_0_0'}, axis=1, inplace=True)
		return df
		
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

	def drop(self, labels, *args, **kwargs):
		"""This not updated for MultiIndex. May not be needed anymore?"""
		labels = self.field2cols(labels)
		return super().drop(labels, *args, **kwargs)

	def field2cols(self, fields):
		"""Find all column names containing 'field' using regex."""
		"""This not updated for MultiIndex. May not be needed anymore?"""
		if isinstance(fields, str):
			fields = [fields]
		fields = [f"f{field}_" for field in fields]
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
		"""This not updated for MultiIndex. May not be needed anymore?"""
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
		"""This not updated for MultiIndex. May not be needed anymore?"""
		cols1 = self.field2cols(field)
		cols2 = self.field2cols(other)
		mask = self[cols2].pkisin(values)
		mask = mask.rename(columns=dict(zip(cols2,cols1)))
		out = self._obj[cols1].mask(~mask, other=pd.NA)
		out = out.mask(out.isin([-1,-3,'-1','-3']), other=pd.NA).dropna(axis='columns', how='all')
		logger.debug(f"findinterpolated: field={field}; other={other}; values={values}; return={out.to_dict()}")
		return out

#	def getfield(self, fields):
#		"""Lookup """
#		self

	def to_psam(self):
		from pkpheno.pkpheno import Psam as Cons
		return self._to_something_multiindex(Cons)

	def to_rvtest(self):
		from pkpheno.pkpheno import RVtest as Cons
		return self._to_something_multiindex(Cons)

	def to_snptest(self, covariates=[]):
		from pkpheno.pkpheno import Snptest as Cons
		return self._to_something_multiindex(Cons)

	def to_textfile(self):
		from pkpheno.textfile import TextFile as Cons
		return self._to_something_multiindex(Cons)

	def _to_something_multiindex(self, cons):
		"""Convert Phenotype Classes must be Overloaded for MultiIndex."""
		obj = self._obj.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None)], level=0, errors='ignore')
		obj.columns = ["_".join(c) for c in obj.columns.to_flat_index()]
		obj[cons.mkey_id] = self.index
		obj[cons.mkey_sex] = self.sex
		try: obj[cons.mkey_altid] = self._obj[self.mkey_altid]
		except AttributeError:
			pass
		textfile = cons(obj)
		return textfile



