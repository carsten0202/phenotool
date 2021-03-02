
#
# ---%%%  pkpheno.py: Handling Phenotype information  %%%---
#

##################################################
#
# --%%	RUN: Perform Basic Setup  %%--

import logging
import numpy as np
import pandas as pd
import sys
import warnings

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import pkcsv as csv

# --%%	END: Perform Basic Setup  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS pheno  %%--

#@pd.api.extensions.register_dataframe_accessor("pheno")
class Phenotype:
	"""A handy-dandy object for storing phenotypes. Will be cool to also import into other scripts when working with phenotypes"""
	__name__ = "Phenotype"
	MAGIC_COLS = {"IID":     ["ID","id","ID_1","Id_1","id_1","IID","Iid","iid","ParticID","Particid","particid"],
	              "FID":     ["FID","Fid","fid","ID_2","Id_2","id_2"],
	              "MAT":     ["MAT"],
	              "MISSING": ["MISSING","Missing","missing"],
	              "PAT":     ["PAT"],
	              "SEX":     ["SEX","Sex","sex","GENDER","Gender","gender"], # SEX should be encoded with males as '1' or 'M'; females as '2' or 'F'
	              "SID":     ["SID"]}

	mkey_id    = "IID" # Also the index, so must be unique.
	mkey_sex   = "SEX"

	def __init__(self, *args, columns=[], samples=[], **kwargs):
		self._obj = pd.DataFrame(*args, **kwargs)
		logger.info(f"{self.__name__}: Parsing file with columns: {self._obj.columns.values}")
		self._obj = self._set_magic_kcol()
		self._obj = self._obj.set_index(self.mkey_id)
		self = self.set_columns(columns=columns)
		self = self.set_samples(samples=samples)
		if self.is_sample_data():
			self._sampletypes = self._obj.index[0]
			self._obj = self._obj.drop(self._obj.index[0])
		Phenotype._validate(self)

# Maybe move this to a method? Like 'standardize', or 'conform'?
		self._obj.replace('NA', np.NaN, inplace=True)
		for col in self._obj.columns:
			self._obj[col] = pd.to_numeric(self._obj[col], errors='ignore')
		self._obj = self._obj.convert_dtypes()

	def _set_magic_kcol(self):
		for k, v in self.MAGIC_COLS.items():
			self._obj = self._obj.rename(dict(zip(v, [k] * len(v))), axis=1)
		assert self.mkey_id in self._obj, f"Unable to identify primary identifier in one or more input files. Please make sure that each input file has a column with one of these recognized headings: {', '.join(self.MAGIC_COLS[self.mkey_id])}"
		return self._obj

	@staticmethod
	def _validate(self, warn=False):
		# Check What? 
		# We should check the SEX column, if it is there... And we need to estimate 0,1,2 errors (when 1 cn mean both males (psam) and females (rvtest)).
		# We should check that all columns have headers and that the headers are strings (numerical headers usually indicates a problem).
		# We should check that there is a header (difficult...)
		# Need to recognize R's one-header-name-short-stupid format
		# Option to set the ID and a warning if it's not recognized
		#	We actually need ways to specify custom names for all magic cols...
		dup = self._obj.index.duplicated()
		assert not dup.any(), f"One or more input files had duplicated primary identifiers ({self._obj.index[dup].unique().to_list()}). The primary identifiers must be unique within each input file."
		if warn:
			# Check for suspect values like 'NA' (Check for Inf?)
			for allna in [s[1].isna().all() for s in self._obj.iterrows()]:
				if allna:
					# SHIT! These guy don't trgger properly if we have columns with auto-default values, like '0'...
					warnings.warn("Dataset contains samples without any associated data; one or more rows will have only missing values.")
					# We can add logging.warning('message') here for additional details.
			for anyna, allna in [(s[1].isna().any(), s[1].isna().all()) for s in self._obj.iteritems()]:
				if allna:
					warnings.warn("Dataset contains phenotypes without any associated data; one or more columns will have only missing values.")
				elif anyna:
					logger.debug("The dataset contains one or more missing values.")

	def is_sample_data(self):
		"""0 (for the first identifier column)
		D (for a column containing discrete values, e.g. a set of strings)
		P or C - for columns containing continuous value - each value must be numerical, or a missing value.
		B - for a column containing a binary trait. The values in this column must be '0', '1', 'control', or 'case'.
		"""
		list_ = [0,"0","B","C","D","P"]
		return all(self._obj.iloc[0].isin(list_))

	def combine_first(self, other):
		"""Return: Combined DataFrame from self and other."""
		if isinstance(other, Phenotype):
			self._obj = self._obj.combine_first(other=other._obj)
		else:
			self._obj = self._obj.combine_first(other=other)
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
		"""Return: All magic columns in input + specified columns in that order;
		columns: a list of column names to extract. Default: All.
		"""
		columns = self.colnames_translate(columns) if columns else self.colnames
		index = pd.Index(self.colnames_magic).append(pd.Index(columns)).drop_duplicates().to_list()
		cols = [col for col in index if col in self._obj]
		self._obj = self._obj.loc[:,cols]
		return self

	def set_samples(self, samples=[]):
		"""Return: All specified samples in specified order."""
		if samples:
			self._obj = self._obj.reindex(samples)
		return self

	def write(self, dest=sys.stdout, *args, **kwargs):
		"""Output self._obj to dest."""
		self._validate(self, warn=True)
		self._obj.to_csv(dest, *args, **kwargs)

# --%%	END: Define CLASS pheno  %%--
#
##################################################




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
	
	https://www.well.ox.ac.uk/~gav/qctool_v2/documentation/sample_file_formats.html
	https://jmarchini.org/file-formats/
	"""
	__name__ = "Snptest"
	MAGIC_COLS = {'ID_1'    : Phenotype.MAGIC_COLS['IID'],
	              'ID_2'    : Phenotype.MAGIC_COLS['FID'],
	              'missing' : Phenotype.MAGIC_COLS['MISSING'],
	              'sex'     : Phenotype.MAGIC_COLS['SEX']}

	mkey_id    = "ID_1" # Also the index, so must be unique.
	mkey_altid = "ID_2"
	mkey_sex   = "sex"

	def __init__(self, *args, columns=[], covariates=[], phenotypes=[], **kwargs):
		""""""
		super().__init__(*args, columns=columns+covariates+phenotypes, **kwargs)
		self._obj[self.mkey_altid] = self._obj.index
		self._obj['missing'] = self._obj.get('missing', pd.Series(pd.NA * self._obj.index.size, index=self._obj.index))
		self._obj[self.mkey_sex] = pd.Categorical(self._obj.get(self.mkey_sex, [""]*len(self._obj.index)))
		self._obj[self.mkey_sex] = self._obj[self.mkey_sex].cat.rename_categories({0: pd.NA, 1: 'male', 2: 'female'})
		self.covariates = covariates
		self.phenotypes = phenotypes
		Snptest._validate(self)

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
# coltype can find them, but they are not getting set as categorical from the stert. Should they be?
# TODO: Use the types from the arguments and input files.
		stype = pd.Series(index=self._obj.columns, dtype='object')
		for colname in self._obj.columns:
			Dtype = self._obj[colname].dropna().convert_dtypes().dtype
			if colname in [self.mkey_id, self.mkey_altid, 'missing']:
				stype[colname] = '0'
			elif Dtype.name == 'category':
				stype[colname] = 'B' if all(self._obj[colname].dropna().isin([0,1,'0','1','Case','case','Control','control'])) else 'D'
			elif pd.api.types.is_numeric_dtype(Dtype):
				stype[colname] = 'C'
			elif Dtype.name == 'string':
				stype[colname] = 'D'
			else:
				warnings.warn("Unable to properly identify datatype for column '{0}'. Setting to 'D' (Discrete values).".format(colname))
		stype[self.covariates] = 'C'
		stype[self.phenotypes] = 'P'
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
		"""Super(), then set covariates/phenotypes."""
		obj = super().combine_first(other, *args, **kwargs)
		obj.covariates = self.covariates.union(other.covariates)
		obj.phenotypes = self.phenotypes.union(other.phenotypes)
		return obj

	def write(self, *args, dest=sys.stdout, **kwargs):
# NOTE: All phenotypes should appear after the covariates in this file.
		self = self.set_columns()
		print(' '.join([self.mkey_id] + self._obj.columns.to_list()))
		print(' '.join(['0'] + self.coltype.to_list()))
		super().write(dest, *args, sep=' ', na_rep='NA', header=False, **kwargs)

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
		""""""
		super().__init__(*args, **kwargs)
		self._obj[self.mkey_id]    = self._obj.index
		self._obj[self.mkey_altid] = self._obj.index
		self._obj[self.mkey_pat]   = self._obj.get(self.mkey_pat, 0)
		self._obj[self.mkey_mat]   = self._obj.get(self.mkey_mat, 0)
		self._obj[self.mkey_sex]   = pd.Categorical(self._obj.get(self.mkey_sex, ['']*len(self._obj.index)))
		self._obj[self.mkey_sex]   = self._obj[self.mkey_sex].cat.rename_categories({0 : pd.NA, 'M' : 1, 'Male' : 1, 'male' : 1, 'F' : 2, 'Female' : 2, 'female' : 2})
		Psam._validate(self)

	def write(self, *args, dest=sys.stdout, header=True, **kwargs):
		"""Like super() but handles *.fam files through header=False if self._obj only has one phenotype."""
		self = self.set_columns()
		if header:
			print("#", end="")
		super().write(dest, *args, sep='\t', na_rep='NA', header=header, index=False, **kwargs)

# --%%	END: Define CLASS Psam (Plink2)  %%--
#
##################################################




##################################################
#
# --%%	RUN: Define CLASS RVtest  %%--

class RVtest(Psam):
	"""Holds RVtest phenotype data
	Documentation can be found here:
	http://zhanxw.github.io/rvtests/#phenotype-file
	"""
	__name__ = "RVtest"
	MAGIC_COLS = {
		"fid":     Phenotype.MAGIC_COLS["FID"],
		"iid":     Phenotype.MAGIC_COLS["IID"], # In Actuality, the ID 'column' isn't a column, it's the index.
		"fatid":   Phenotype.MAGIC_COLS["MAT"],
		"matid":   Phenotype.MAGIC_COLS["PAT"],
		"sex":     Phenotype.MAGIC_COLS['SEX'], # SEX should be encoded with males as '0', females as '1' and missing as 'NA'
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

	def write(self, *args, dest=sys.stdout, **kwargs):
		"""Super() skips over Psam to avoid the initial '#'."""
		self = self.set_columns()
		super(Psam,RVtest).write(self, dest, *args, sep='\t', na_rep='NA', index=False, **kwargs)

# --%%	END: Define CLASS RVtest  %%--
#
##################################################




##################################################
#
# --%%  RUN: Define CLASS PEP for Portable Encapsulated Projects  %%--

class PEP(Phenotype):
	"""Portable Encapsulated Projects

	"""
	__name__ = "PEP"

	def __init__(self, *args, **kwargs):
		""""""
		super().__init__(*args, **kwargs)
		PEP._validate(self)

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
	Inspired by
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

