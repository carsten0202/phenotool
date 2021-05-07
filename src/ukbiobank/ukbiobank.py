

import logging
#import numpy as np
import pandas as pd
import re
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

#import pklib.pkcsv as csv
from pkpheno import Phenotype, Psam



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
		"SEX":     ["f.31.0.0","31-0.0"], # UKB encodes as 0=female, but SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
	}

	mkey_id    = "EID" # Also the index, so must be unique.
	mkey_sex   = "SEX"

	def __init__(self, iterable, *args, phenovars, samples=[], **kwargs):
		"""
		iterable: An iterable with data...
		phenovars: The UKBiobank datafield(s) to extract. (Named for compatibility with ancestor classes).
		"""
# NOTE: The second digit in datafields is called an 'instance'.
# NOTE: The third is the array index.
		def genfunc(dicttable, cols, rows=[], indexcol=0):
			"""Will extract cols from a iterable dict table. Cols will be outputted only once regardless of duplicates"""
			for dictrow in dicttable:
				row = list(dictrow.values())
				if not rows or row[indexcol] in rows:
					yield [row[col] for col in list(dict.fromkeys(cols))]

		icol = [0]
		headdict = next(iterable)
		for field in phenovars + [31]:
			found = False
			for i, col in enumerate(headdict):
				if re.match(f"\D*{field}\D", col) is not None:
					icol.append(i)
					found = True
			if not found:
				logger.error(f"UKBioBank: Datafield '{field}' does not exist in input data.")
		logger.debug(f"UKBioBank: Selected columns {[list(headdict.keys())[i] for i in list(dict.fromkeys(icol))]}")
		super().__init__(genfunc(iterable, icol, rows=samples), *args, columns=[list(headdict.keys())[i] for i in list(dict.fromkeys(icol))], samples=samples, **kwargs)
		patone = re.compile("f\.")
		pattwo = re.compile("[-.]")
		self._obj = self._obj.rename(columns=lambda label: pattwo.sub("_",patone.sub('f',label)))
		logger.debug(f"UKBioBank: Renamed columns {list(self._obj.columns)}")

	@property
	def sex(self):
		"""Returns the SEX in a systematic way (male/female) for querying."""
		out = pd.Series([pd.NA] * self._obj.index.size, name=self.mkey_sex, index=self._obj.index)
		out[self[self.mkey_sex].pkisin(['1'])] = "male"
		out[self[self.mkey_sex].pkisin(['0'])] = "female"
		return out

	def drop(self, labels, *args, **kwargs):
		labels = self.field2cols(labels)
		return super().drop(labels, *args, **kwargs)

	def findinfield(self, fields, values):
		"""Find value in indicated fields.

		Return:
		    A pd.Series with True/False for each row describing if row contained the value.

		Note: It seems -1 and -3 are consistently used to indicate 'missing' in UKB. This is implemented here.
		"""
		# Ok, here's a serious bug. Value 'NA' in UKB doesn't mean 'pd.NA', rather it means False.
		if isinstance(fields, str):
			fields = [fields]
		fields = [f"f{field}_" for field in fields]
		logger.debug(f"findinfield: Scanning fields='{fields}' for values={values}.")
		out = pd.Series(pd.NA, index=self.index, dtype=pd.BooleanDtype())
		cols = self.field2cols(fields)
		if self.columns.isin(cols).any():
			out[self.index] = self[cols].pkisin(values).any(axis='columns')
			out[self._obj[cols].mask(self[cols].pkisin(['-1','-3']), pd.NA).isna().all(axis='columns')] = pd.NA
		logger.debug(f"findinfield: Found = {out.to_dict()}.")
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

	def findfield(self, fields):
		"""Return any and all values in a given field."""
		return self._obj[self.field2cols(fields)]

