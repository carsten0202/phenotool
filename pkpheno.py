#!/home/fls530/anaconda3/bin/python

# Module to hold the Pheno class and its associated functions

##################################################
#
# --%%	RUN: Perform Basic Setup  %%--

import numpy as np
import pandas as pd
import pkcsv as csv
import sys
import warnings

# --%%	END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS pheno  %%--

#@pd.api.extensions.register_dataframe_accessor("pheno")
class Phenotype:
	"""A handy-dandy object for storing phenotypes. Will be cool to also import into other scripts when working with phenotypes"""

	MAGIC_COLS = {"IID":     ["ID","id","ID_1","Id_1","id_1","IID","Iid","iid"], # Not really magic; the ID column is always the first column, regardless of name
	              "FID":     ["FID","Fid","fid","ID_2","Id_2","id_2"],
	              "MAT":     ["MAT"],
	              "MISSING": ["MISSING","Missing","missing"],
	              "PAT":     ["PAT"],
	              "SEX":     ["SEX","Sex","sex"], # SEX should be encoded with males as '1' or 'M'; females as '2' or 'F'
	              "SID":     ["SID"]}

	def __init__(self, *args, **kwargs):
		self._obj = pd.DataFrame(*args, **kwargs)
		self._validate(self._obj)
		self._obj = self._set_magic_kcol()
		self._obj.replace('NA', np.NaN, inplace=True)
		binaries = []
		if self.is_sample_data():
			binaries = self._obj.columns[['B' == v for v in self._obj.iloc[0].values]]
			self._obj = self._obj.drop(self._obj.index[0])
		self._obj = self._obj.set_index(self._obj.columns[0])
		for col in self._obj.columns:
			self._obj[col] = pd.to_numeric(self._obj[col], errors='ignore')
		self._obj = self._obj.convert_dtypes()
		self._obj['SEX'] = pd.Categorical(self._obj['SEX'])
		for cat in binaries:
			self._obj[cat] = pd.Categorical(self._obj[cat])

	def _set_magic_kcol(self):
		for k, v in self.MAGIC_COLS.items():
			self._obj = self._obj.rename(dict(zip(v, [k] * len(v))), axis=1)
		return self._obj

	@staticmethod
	def _validate(obj):
		# Check What? 
		# We should check the SEX column, if it is there...
		# We should check that there is a header (difficult...)
		# We should check that there is an id column
		# We could check for 'Inf' and other suspect values
		# We could check that there are no duplicate magic cols? No duplicates in general?
		# Sex? ()
		pass

	def is_sample_data(self):
		"""0 (for the first identifier column)
		D (for a column containing discrete values, e.g. a set of strings)
		P or C - for columns containing continuous value - each value must be numerical, or a missing value.
		B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'."""
		list_ = [0,"0","B","C","D","P"]
		return all(self._obj.iloc[0].isin(list_))

	def combine_first(self, other):
		"""Return: Combined DataFrame from self and other. NOTE: Will fail if self and other are not both Phenoypes!"""
		self._obj = self._obj.combine_first(other=other._obj)
		self._obj["SEX"] = pd.Categorical(self._obj["SEX"])
		return self

	def get_all_magic(self):
		"""Return: All magic cols in self as list"""
		mcols = [self.get_magic_kcol(mcol) for mcol in self.MAGIC_COLS if self.get_magic_kcol(mcol)]
		return pd.Index(mcols).drop_duplicates().to_list()

	def get_magic_kcol(self, colname):
		"""Return: Names of magic cols as list or str (same type as 'colname')"""
		if isinstance(colname, list):
			mcols = []
			for col in colname:
				mcols.extend([k for k in self.MAGIC_COLS[col] if k in self._obj][0:1])
			return mcols
		else:
			for kcol in [k for k in self.MAGIC_COLS[colname] if k in self._obj]:
				return kcol
		return None

	def get_phenotypes(self, phenotypes=None):
		"""Return: Specified columns + all magic columns; columns: a list of column names to extract"""
# This will fail if phenotypes contains magic cols under alternative names

		index = pd.Index(self.get_all_magic() + phenotypes).drop_duplicates().to_list()
		cols = [col for col in index if col in self._obj]
		self._obj = self._obj.loc[:,cols]
		if not index in cols:
			warnings.warn("Not all specified phenotypes were found in input files.")
		return self

	def get_samples(self, samples):
# Arghhh.. it fails bedause all samples are numeric??!!!
		self._obj = self._obj.reindex(samples)
		# INSERT: warnings.warn("Not all samples were found in input. Some subjects will have all missing values.")
		return self

# --%%	END: Define CLASS pheno  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS Snptest  %%--

#@pd.api.extensions.register_dataframe_accessor("pheno")
class Snptest(Phenotype):
	"""Holds phenotypes in a format appropriate for Snptest"""

	MAGIC_COLS = {'ID_2'    : Phenotype.MAGIC_COLS['FID'],
	              'missing' : Phenotype.MAGIC_COLS['MISSING'],
	              'sex'     : Phenotype.MAGIC_COLS['SEX']}

	def __init__(self, obj, *args, **kwargs):
		if isinstance(obj, Phenotype):
			self._obj = obj._obj
			self._obj = self._set_magic_kcol()
		else:
			super().__init__(obj, *args, **kwargs)
		# We should then do some coversion here..
		self._obj['ID_2'] = self._obj.index
		self._obj["sex"] = self._obj["sex"].cat.rename_categories({1: 'MALE', 2: 'FEMALE'})

#	@staticmethod
#	def _validate(obj):
		# Check What? 
		# We should check that there is a header (difficult...)
		# We should check that there is an id column
		# We could check for 'Inf' and other suspect values
		# We could check that there are no duplicate magic cols? No duplicates in general?
		# Sex? ()
#		super()._validate(obj)

	def get_sample_type(self):
		"""Set the sample type for the pseudo-type-header in sample files.
		Values are:
		0 (for the first identifier column)
		D (for a column containing discrete values, e.g. a set of strings)
		P or C - for columns containing continuous value - each value must be numerical, or a missing value.
		B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.
		Sex of samples should be in column 'sex' of type D. Values in the column must be "f", "m", "male" or "female".
		"""
# WAANING: No proper support for binary traits yet... should be 'control', or 'case'.
		stype = pd.Series(index=self._obj.columns, dtype='object')
		for colname in self._obj.columns:
			Dtype = self._obj[colname].dropna().convert_dtypes().dtype
			if Dtype.name == 'category':
				stype[colname] = 'B' if len(self._obj[colname].cat.categories) <= 2 else 'D'
			elif pd.api.types.is_numeric_dtype(Dtype):
				stype[colname] = 'C'
			elif Dtype == 'string':
				stype[colname] = 'D'
			else:
				warnings.warn("Unable to properly identify datatype for {0}. Setting to 'D'.".format(colname))
		stype['ID_2'] = '0'
		stype['missing'] = '0'
		return stype

	def write(self):
		neworder = pd.Index(["ID_2","missing","sex"]).append(self._obj.columns).drop_duplicates(keep='first')
		self._obj = self._obj[neworder]
		print(' '.join(['ID_1'] + self._obj.columns.to_list()))
		print(' '.join(['0'] + self.get_sample_type().to_list()))
		self._obj.to_csv(sys.stdout, sep=' ', na_rep='NA', header=False)

# --%%	END: Define CLASS Snptest  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS Psam (Plink2)  %%--

class Psam(Phenotype):
	"""Holds phenotypes in a format appropriate for Plink2"""
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	def writepsam(self):
		self._obj["IID"] = self._obj.index
		if self._get_magic_kcol("MISSING"):
			self._obj.drop(columns=self._get_magic_kcol("MISSING"), inplace=True)
		mcols = {self._get_magic_kcol(c):c for c in ["FID","IID","SID","PAT","MAT","SEX"] if self._get_magic_kcol(c)}
		self._obj.rename(columns=mcols,inplace=True)
		neworder = pd.Index(mcols.values()).append(self._obj.columns).drop_duplicates(keep='first')
		self._obj = self._obj[neworder]
		print("#", end="")
		self._obj.to_csv(sys.stdout, sep='\t', na_rep='NA', index=False)

# --%%	END: Define CLASS Psam (Plink2)  %%--
#
##################################################




##################################################
#
# --%%  RUN: Constructors  %%--

def fromiter(phenoiter, *args, phenotypes=None, samples=None, **kwargs):
	df = Phenotype(list(phenoiter), *args, **kwargs)
	if phenotypes is not None:
		df = df.get_phenotypes(phenotypes=phenotypes)
	if samples is not None:
		df = df.get_samples(samples=samples)
	return df

# --%% END: Constructors  %%--
#
##################################################

