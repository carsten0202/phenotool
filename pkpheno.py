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

	def __init__(self, *args, columns=[], samples=[], **kwargs):
		self._obj = pd.DataFrame(*args, **kwargs)
		self._obj = self._obj.set_index(self._obj.columns[0])
		self._obj = self._set_magic_kcol()
		self._validate(self)
		self = self.set_columns(columns=columns)
		self = self.set_samples(samples=samples)
		if self.is_sample_data():
			self._sampletypes = self._obj.index[0]
			self._obj = self._obj.drop(self._obj.index[0])

# Maybe move this to a method? Like 'standardize', or 'conform'?
		self._obj.replace('NA', np.NaN, inplace=True)
		for col in self._obj.columns:
			self._obj[col] = pd.to_numeric(self._obj[col], errors='ignore')
		self._obj = self._obj.convert_dtypes()

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
		self._obj = self._obj.fillna(np.NaN)
		return self

	@property
	def colnames(self):
		"""Return: All colnames in _obj as list"""
		return self._obj.columns.to_list()

	@property
	def colnames_magic(self):
		"""Return: All magic colnames supported by class as list."""
		return list(self.MAGIC_COLS.keys())

	@property
	def colnames_normal(self):
		"""Return: All colnames in _obj which are not magical."""
		return [col for col in self._obj.columns if col not in self.colnames_magic]

	def columns_rankINT(self, columns=[]):
		"""Return: Self after Rank-based inverse normally transformation on cols given by 'columns'. Default: All normal columns"""
		columns = self.colnames_translate(columns) if columns else self.colnames_normal
		for col in columns:
			self._obj[col] = rank_INT(self._obj[col])
		return self

	def colnames_translate(self, columns):
		"""Translate names in columns to the canonical magic column names used by self."""
		outcols = []
		for col in columns:
			mcol = [k for k in self.colnames_magic if col in self.MAGIC_COLS[k]]
			outcols.extend(mcol if mcol else [col])
		return outcols

	def set_columns(self, columns=[]):
		"""Return: Specified columns + all magic columns; columns: a list of column names to extract"""
		columns = self.colnames_translate(columns) if columns else self.colnames
		index = pd.Index(self.colnames_magic).append(pd.Index(columns)).drop_duplicates().to_list()
		cols = [col for col in index if col in self._obj]
		self._obj = self._obj.loc[:,cols]
		if not index in cols:
			# This warning triggers too often/early
			warnings.warn("Not all specified phenotypes were found in input files.")
		return self

	def set_samples(self, samples=[]):
		if samples:
			self._obj = self._obj.reindex(samples)
		# INSERT: warnings.warn("Not all samples were found in input. Some subjects will have all missing values.")
		return self

	def write(self):
		print(' '.join(['ID'] + self.colnames))
		self._obj.to_csv(sys.stdout, sep=' ', na_rep='NA', header=False)

# --%%	END: Define CLASS pheno  %%--
#
##################################################



# Inspired by: https://github.com/edm1/rank-based-INT/blob/master/rank_based_inverse_normal_transformation.py
def rank_INT(series, c=3.0/8, stochastic=True):
	"""Perform rank-based inverse normal transformation on pandas series.
	   If stochastic is True ties are given rank randomly, otherwise ties will
	   share the same value. NaN values are ignored.
	   Args:
		param1 (pandas.Series):   Series of values to transform
		param2 (Optional[float]): Constand parameter (Bloms constant)
		param3 (Optional[bool]):  Whether to randomise rank of ties
	   Returns:
		pandas.Series"""
	import os
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



##################################################
#
# --%%	RUN: Define CLASS Snptest  %%--

# TODO: We currently have no way to properly calculate the missing data proportion of each individual. (The 'missing' Column)
# TODO: Ok, we need to specify in detail the phenotypes and covariates. 
#	Put them in __init__? self.covar = stuff ?
#	Also means we must read files directly as Snptest....

#@pd.api.extensions.register_dataframe_accessor("pheno")
class Snptest(Phenotype):
	"""Holds phenotypes in a format appropriate for Snptest
	Documentation can be found here:
	
	https://jmarchini.org/file-formats/
	"""

	MAGIC_COLS = {'ID_2'    : Phenotype.MAGIC_COLS['FID'],
	              'missing' : Phenotype.MAGIC_COLS['MISSING'],
	              'sex'     : Phenotype.MAGIC_COLS['SEX']}

	def __init__(self, *args, columns=[], covariates=[], phenotypes=[], **kwargs):
		super().__init__(*args, columns=columns+covariates+phenotypes,**kwargs)
		# We should then do some coversion here..
		self._obj['ID_2'] = self._obj.index
		self._obj['sex'] = pd.Categorical(self._obj.get('sex', [""]*len(self._obj.index)))
		self._obj["sex"] = self._obj["sex"].cat.rename_categories({1: 'male', 2: 'female'})
		self.covariates = covariates
		self.phenotypes = phenotypes

	@property
	def coltype(self):
		"""Set the column type for the pseudo-type-header in sample files.
		Values are:
		0 (for the first identifier column)
		D (for a column containing discrete values, e.g. a set of strings)
		P or C - for columns containing continuous value - each value must be numerical, or a missing value.
		B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.
		Sex of samples should be in column 'sex' of type D. Values in the column must be "f", "m", "male" or "female".
		"""
# WARNING: No proper support for binary traits yet... must only be '0 = control', or '1 = case'
# coltype can find them, but they are not getting set as categorical from the stert. Should they be?
# TODO: Use the types from the arguments and input files.
		stype = pd.Series(index=self._obj.columns, dtype='object')
		for colname in self._obj.columns:
			Dtype = self._obj[colname].dropna().convert_dtypes().dtype
			if Dtype.name == 'category':
				stype[colname] = 'B' if all(self._obj[colname].dropna().isin([0,1,'0','1','case','control'])) else 'D'
			elif pd.api.types.is_numeric_dtype(Dtype):
				stype[colname] = 'C'
			elif Dtype == 'string':
				stype[colname] = 'D'
			else:
				warnings.warn("Unable to properly identify datatype for {0}. Setting to 'D'.".format(colname))
		stype[self.covariates] = 'C'
		stype[self.phenotypes] = 'P'
		stype['ID_2'] = '0'
		stype['missing'] = '0'
		return stype

	@property
	def covariates(self):
		return self._covariates

	@property
	def phenotypes(self):
		return self._phenotypes

	@covariates.setter
	def covariates(self, value):
		value = pd.Index(value)
		self._covariates = value.intersection(self._obj.select_dtypes(include='number').columns)

	@phenotypes.setter
	def phenotypes(self, value):
		value = pd.Index(value)
		self._phenotypes = value.intersection(self._obj.select_dtypes(include='number').columns)

	def combine_first(self, other, *args, **kwargs):
		obj = super().combine_first(other, *args, **kwargs)
		obj.covariates = self.covariates.union(other.covariates)
		obj.phenotypes = self.phenotypes.union(other.phenotypes)
		return obj

	def write(self):
# NOTE: All phenotypes should appear after the covariates in this file.
		neworder = pd.Index(["ID_2","missing","sex"]).append(self._obj.columns).drop_duplicates(keep='first')
		self._obj = self._obj[neworder]
		print(' '.join(['ID_1'] + self._obj.columns.to_list()))
		print(' '.join(['0'] + self.coltype.to_list()))
		self._obj.to_csv(sys.stdout, sep=' ', na_rep='NA', header=False)

# --%%	END: Define CLASS Snptest  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS Psam (Plink2)  %%--

class Psam(Phenotype):
	"""Holds phenotypes in a format appropriate for Plink2.
	Documentation can be found here:
	https://www.cog-genomics.org/plink/2.0/formats#psam
	"""

	MAGIC_COLS = {"IID":     Phenotype.MAGIC_COLS["IID"], # Not really magic; the ID column is always the first column, regardless of name
	              "FID":     Phenotype.MAGIC_COLS["FID"],
	              "MAT":     Phenotype.MAGIC_COLS["MAT"],
	              "PAT":     Phenotype.MAGIC_COLS["PAT"],
	              "SEX":     Phenotype.MAGIC_COLS['SEX'], # SEX should be encoded with males as '1', females as '2' and missing as '0' or 'NA'
	              "SID":     Phenotype.MAGIC_COLS['SID']} # SID is not really supported

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		self._obj['IID'] = self._obj.index
		self._obj['PAT'] = self._obj.get('PAT', 0)
		self._obj['MAT'] = self._obj.get('MAT', 0)
		self._obj['SEX'] = pd.Categorical(self._obj.get('SEX', ['']*len(self._obj.index)))
		self._obj['SEX'] = self._obj['SEX'].cat.rename_categories({'M' : 1, 'Male' : 1, 'male' : 1, 'F' : 2, 'Female' : 2, 'female' : 2})

	def combine_first(self, other, *args, **kwargs):
		obj = super().combine_first(other, *args, **kwargs)
		return obj

	def write(self):
		neworder = pd.Index(["FID","IID","PAT","MAT","SEX"]).append(self._obj.columns).drop_duplicates(keep='first')
		self._obj = self._obj[neworder]
		print("#", end="")
		self._obj.to_csv(sys.stdout, sep='\t', na_rep='NA', index=False)

# --%%	END: Define CLASS Psam (Plink2)  %%--
#
##################################################




##################################################
#
# --%%  RUN: Constructors  %%--


# --%% END: Constructors  %%--
#
##################################################

