
###########################################################
#
# ---%%%  pkpheno.py: Handling Phenotype information  %%%---
#

##################################################
#
# --%%	RUN: Perform Basic Setup  %%--

import copy
import logging
from numbers import Number
import numpy as np
import pandas as pd
import sys
import warnings

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

import pklib.pkcsv as csv

# --%%	END: Perform Basic Setup  %%--
#
##################################################



##################################################
#
# --%%	RUN: Define CLASS pheno  %%--

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

	def __init__(self, *args, phenovars=[], samples=[], **kwargs):
		self._obj = pd.DataFrame(*args, **kwargs)
		logger.info(f"{self.__name__}: Parsing file with columns[:10] = {self._obj.columns.to_list()[:10]}")
		logger.debug(f"{self.__name__}: Parsing file with columns: {self._obj.columns.to_list()}")
		self._obj = self._set_magic_kcol()
		self._obj = self._obj.set_index(self.mkey_id)
		self = self.set_columns(columns=phenovars)
		if samples:
			self.samples = samples
		if self.is_sample_data():
			self._sampletypes = self._obj.index[0]
			self._obj = self._obj.drop(self._obj.index[0])
		Phenotype._validate(self)

# Maybe move this to a method? Like 'standardize', or 'conform'?
		self._obj.replace('NA', np.NaN, inplace=True)
		for col in self._obj.select_dtypes('object').columns:
			self._obj[col] = pd.to_numeric(self._obj[col], errors='ignore')
		self._obj = self._obj.convert_dtypes()
		logger.debug(f"Columns after __init__: {self.columns.to_list()}")

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
			self._obj = self._obj.rename(dict(zip(v, [k] * len(v))), axis=1)
		assert self.mkey_id in self._obj, f"Unable to identify primary identifier in one or more input files. Please make sure that each input file has a column with one of these recognized headings: {', '.join(self.MAGIC_COLS[self.mkey_id])}"
		return self._obj

	@staticmethod
	def _validate(self, warn=False):
		# Check What? 
		# We should check the SEX column, if it is there... And we need to estimate 0,1,2 errors (when 1 cn mean both males (psam) and females (rvtest)).
		# We should check that all columns have headers and that the headers are strings (numerical headers usually indicates a problem).
		# We should check that there is a header (difficult...)
		dup = self._obj.index.duplicated()
		assert not dup.any(), f"One or more input files had duplicated primary identifiers ({self._obj.index[dup].unique().to_list()}). The primary identifiers must be unique within each input file."
		if warn:
			samples_notok = pd.Series(pd.NA * self.index.size, index=self.index)
			columns_notok = pd.Series(pd.NA * len(self.columns), index=self.columns)
			for name,row in self.df.iterrows():
				if row.isna().all():
					samples_notok[name] = name
			samples_notok = samples_notok.dropna().to_list()
			for name,col in self.df.iteritems():
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
		return [col for col in self._obj.columns if col not in self.colnames_magic]

	@property
	def columns(self):
		"""Returns the columns of _obj."""
		return self._obj.columns

	@property
	def df(self):
		"""Return the pheno DataFrame."""
		return self._obj

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
		try: self._obj = self._obj.reindex(value)
		except Exception as ex:
			sys.exit(f"{ex} in samples.setter...")

	@property
	def sex(self):
		"""Returns the SEX in a systematic way (male/female or child-dependent) for querying."""
		out = pd.Series(np.zeros(self._obj.index.size) + np.nan, name=self.mkey_sex, index=self._obj.index)
		out[self[self.mkey_sex].pkisin(['1', 'M'])] = "male"
		out[self[self.mkey_sex].pkisin(['2', 'F'])] = "female"
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
			self._obj = self._obj.combine_first(other=other._obj)
		else:
			self._obj = self._obj.combine_first(other=other)
		self._obj = self._obj.fillna(np.NaN)
		return self

	def derive_rankINT(self, columns=None):
		"""Return: DataFrame with results from a rank-based inverse normally transformation on cols given by
		'columns'.
		Default: All normal columns."""
		columns = self.colnames_translate(columns) if columns else self.colnames_normal
		df = pd.DataFrame()
		for col in columns:
			df[col] = rank_INT(self._obj[col])
		return df

	def drop(self, *args, **kwargs):
		"""NOTE: This guy actually drops inplace..."""
		self._obj = self._obj.drop(*args, **kwargs)
		return self

	def field2cols(self, fields):
		"""Find all column names containing 'field' using regex."""
		out = pd.Series()
		cols = pd.Series(self._obj.columns)
		if not isinstance(fields, list):
			fields = [fields]
		for field in fields:
			out = out.append(cols[cols.str.contains(str(field), regex=True)])
		logger.debug(f"field2cols: Fields={fields}; Found cols={out.to_list()}.")
		return out

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

	def pkisin(self, values):
		"""Convienience function to check for 'values' using both numeric and string in df.isin()."""
		if isinstance(values, str):
			values = [values]
		out = values
		for value in values:
			if isinstance(value, str) and value.isnumeric():
				out.append(float(value))
			elif isinstance(value, Number):
				out.append(str(value))
		return self._obj.isin(out)

	def set_columns(self, columns=[]):
		"""Return: All magic columns in input + specified columns in that order;
		columns: a list of column names to extract. Default: All.
		"""
		columns = self.colnames_translate(columns) if columns else self.colnames
		index = pd.Index(self.colnames_magic).append(pd.Index(columns)).drop_duplicates().to_list()
		cols = [col for col in index if col in self._obj]
		self._obj = self._obj.loc[:,cols]
		return self

	def to_psam(self):
		"""Convert Phenotype Class to Class Psam for Plink output."""
		obj = self._obj.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None), getattr(self, 'mkey_mat', None), getattr(self, 'mkey_pat', None)], errors='ignore')
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
		obj = self._obj.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None)], errors='ignore')
		obj[Snptest.mkey_id] = self.index
		obj[Snptest.mkey_sex] = self.sex
		try: obj[Snptest.mkey_altid] = self._obj[self.mkey_altid]
		except AttributeError:
			pass
		snptest = Snptest(obj, covariates=covariates)
		return snptest

	def to_textfile(self):
		"""Convert Phenotype Class to Class TextFile for custom file output."""
		from .textfile import TextFile
		obj = self._obj.drop(columns=[self.mkey_id, self.mkey_sex, getattr(self, 'mkey_altid', None)], errors='ignore')
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
		self._covariates = value.intersection(self._obj.select_dtypes(include='number').columns).to_list()

	@property
	def phenotypes(self):
		return self._phenotypes

	@phenotypes.setter
	def phenotypes(self, value):
		value = pd.Index(value)
		self._phenotypes = value.intersection(self._obj.select_dtypes(include='number').columns)

# Probably don't need these guys: the Phenotype default should suffice
#	@property
#	def sex(self):
#		"""Call getter from parent."""
#		return super().sex

#	@sex.setter
#	def sex(self, value):
#		"""Set the sex column."""
#		super(Snptest, self.__class__).sex.fset(self, value)

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

	def to_psam(self):
		"""No conversion necessary; self is already Psam."""
		return self

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

	def to_rvtest(self):
		"""No conversion necessary; self is already RVtest."""
		return self

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
	MAGIC_COLS = {
		"sample_name":     Phenotype.MAGIC_COLS["IID"] + ["sample_name"], # In Actuality, the ID 'column' isn't a column, it's the index.
	}
	mkey_id    = "sample_name" # Also the index, so must be unique.

	def __init__(self, *args, **kwargs):
		""""""
		super().__init__(*args, **kwargs)
		PEP._validate(self)

	def write(self, *args, quoting=csv.QUOTE_ALL, **kwargs):
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

