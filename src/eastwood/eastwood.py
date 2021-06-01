

###########################################################
#
# --%%  Setup and Initialize  %%--

import copy
from datetime import datetime
import logging
import pandas as pd
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)


# Convienience function.

# Need to double check this guy for Kleene logic (It doesn't follow it completely)
def pkall(df, axis=0):
	"""This function is necessary to address bugs with skipna=False when that behavior is needed."""
	out = pd.Series(df.all(axis=axis), dtype='boolean')
	out[df.isna().any(axis=axis)] = pd.NA
	return out

def pkany():
	"""This function is necessary to address bugs with skipna=False when that behavior is needed."""
	sys.exit("Not implemented yet")

def datefirst(df, axis='index'):
	"""Find the earliest occurrance in a pd.Series/pd.DataFrame of dates. Can be applied over a specific axis"""
	out = pd.Series(name="datefirst", dtype='datetime64[ns]', index=df.index)
	for row in df.dropna(axis='index', how='all').itertuples():
		out[row[0]] = pd.to_datetime(row[1:]).min()
	return out




###########################################################
#
# --%%  DEFINE: Eastwood class from the Eastwood2016 paper.  %%--

class Eastwood():
	"""Implements a Diabetes classification based on Eastwood2016."""

	# Define our categories
	CategoriesDM = {
	    'Negative'   : "Diabetes_Unlikely",
	    'T1High'     : "Probable_Type1",
	    'T2High'     : "Probable_Type2",
	    'T1Moderate' : "Possible_Type1",
	    'T2Moderate' : "Possible_Type2",
	    'GDModerate' : "Possible_Gestational",
	}

	# UKB base date
	UKBbaseline = pd.to_datetime("2010-08-01")
	UKBioFields = ['41270', '41280']

	def __init__(self, pheno, baseline=None, datediag=None, **kwargs):
		"""Prevalence and Incidence based on the Eastwood2016 paper."""

		# Start with some assessments. Do we have the data?
		#	Remember that we may want to use this definition on other datasets; not just UKB.
		#		Will make assessment critical, since data could be missing (and user will need to provide much information...)

		# Initialize attributes
		for (k,v) in Eastwood.CategoriesDM.items():
			setattr(self, k, v)
		self.baseline = baseline
		self.dm = pd.DataFrame(index=pheno.index, dtype='boolean')
		self._incidence = pd.Series(
		    [pd.NA] * pheno.index.size,
		    dtype="datetime64[ns]",
		    index=pheno.index,
		    name="Incidence",
		) # This guy is a frozen init; actual incidence is returned through a property getter
		self._prevalence = pd.Series(
		    pd.Categorical([pd.NA] * pheno.index.size,
		                   categories=[self.Negative, self.T1High, self.T2High, self.T1Moderate, self.T2Moderate, self.GDModerate]),
		    index=pheno.index,
		    name="Prevalence",
		) # This guy is a frozen init; actual prevalence is returned through a property getter

		# Needed for both Incidence and Prevalence
		self.dm['date_t1dm_ip']    = datefirst(pheno.findinterpolated('41280', '41270', ['E10' + str(i) for i in range(10)]), axis='columns')
		self.dm['date_t2dm_ip']    = datefirst(pheno.findinterpolated('41280', '41270', ['E11' + str(i) for i in range(10)]), axis='columns')
		self.dm['date_otherdm_ip'] = datefirst(pheno.findinterpolated('41280', '41270', ['E13' + str(i) for i in range(10)]), axis='columns')
		self.dm['date_unkdm_ip']   = datefirst(pheno.findinterpolated('41280', '41270', ['E14' + str(i) for i in range(10)]), axis='columns')
		self.dm['date_anydm_ip']   = datefirst(self.dm[['date_t1dm_ip','date_t2dm_ip', 'date_otherdm_ip', 'date_unkdm_ip']], axis='columns')

		# Validation
		Eastwood._validate(self)

	def __getitem__(self, key):
		"""Propagates index operations across the relevant attributes."""
		out = copy.deepcopy(self)
		out.dm = self.dm.loc[key,]
		out._prevalence = self._prevalence[key]
		return out

	@staticmethod
	def _validate(self):
		"""Check stuff. But what?

		If a parameter matches all or none of the subjects?
		Check that baseline date is valid. For prevalence it should be UKB baseline or later. (what if data are not UKB?)
		"""
		logger.debug(f"validate: dm.columns = {self.dm.columns.to_list()}")
		logger.debug(f"validate: dm.dtypes  = {self.dm.dtypes.to_list()}")

	@property
	def baseline(self):
		"""Get baseline attribute. Mostly here because we need the setter"""
		return self._baseline

	@baseline.setter
	def baseline(self, value):
		"""Baseline setter. Needed for class decendants."""
		self._baseline = value



###########################################################
#
# --%%  DEFINE: Incidence class from the Eastwood2016 paper.  %%--

class Prevalence(Eastwood):
	""""Implements prevalence algorithm 'A+B+C' from Eastwood2016 paper (See Fig 2 in paper).

	The algorithm is slightly modified from the original in that the baseline date is user-adjustable and
	in-patient data prior to baseline are considered.
	"""
	styles = ["Eastwood", 'T1D', 'T2D']
	binaryNegative = 0
	binaryPositive = 1
	binaryOther = pd.NA

	UKBstartdate = pd.to_datetime("2006-01-01")
	UKBioFields = Eastwood.UKBioFields + ['53', '2443', '2976', '2986', '4041', '6153', '6177', '20002', '20003', '20009', '21000']

	DEBUGsubject = '1000863'

	def __init__(self, pheno, *args, agediag=None, ethnicity=None, high=None, moderate=None, style='Eastwood', treatments=None, **kwargs):
		"""Based on the Eastwood2016 paper.

		The following refers to the prevalence algorithm only.
		"""
		super().__init__(pheno, *args, **kwargs)
		self.style = style

		instances = pd.Series(None, index=pheno.index, dtype='object')
		df = pheno.findfield('53').astype('object')
		for (i, s) in (df.astype("datetime64[ns]") <= self.baseline).iterrows():
			instances[i] = df.columns[s].to_series().str.replace(r'f\d+_(\d+)_\d+', r'\1', regex=True).to_list()
		logger.debug(f"Init: Identified instances = {instances.to_list()}")

		# Fill dm with values - Main ethnic groups (NB: this makes some assumptions i.e. British is white, and does not include mixed, Chinese or other Asian)
		self.dm.loc[pheno.findinfield('21000', ['1', '1001', '1002', '1003']), 'ethnic'] = 1 # White European
		self.dm.loc[pheno.findinfield('21000', ['3', '3001', '3002', '3003']), 'ethnic'] = 2 # South Asian
		self.dm.loc[pheno.findinfield('21000', ['4', '4001', '4002', '4003']), 'ethnic'] = 3 # African Caribbean
		self.dm.loc[pheno.findinfield('21000', ['2', '2001', '2002', '2003', '2004', '5', '6']), 'ethnic'] = 4 # Mixed or Other
		self.dm['ethnic_sa_afc'] = pd.Series([ethnic in [2, 3] for ethnic in self.dm['ethnic']], dtype='boolean') # SA and AFC vs all other variable
		logger.info(f"Init: Ethnic breakdown = {dict(zip(['White European', 'South Asian', 'African Caribbean', 'Mixed or Other'], self.dm['ethnic'].value_counts().sort_index()))}.")

		# Fill dm with values - Touchscreen
		self.dm['gdmonly_sr'] = pkall(pd.concat([pheno.findinfield('4041','1',instances), pheno.sex == 'female'], axis='columns'), axis='columns')
		logger.info(f"Init: {self.dm['gdmonly_sr'].sum()} subjects with Gestational Diabetes from Touchscreen.")
#		logger.debug(f"Init: Subject={self.DEBUGsubject} {self.dm.loc[self.DEBUGsubject,]}")

		# Fill dm with values - Nurse Interview
		self.dm['alldm_ni'] = pheno.findinfield('20002','1220',instances)
		logger.info(f"Init: {self.dm['alldm_ni'].sum()} Subjects with any type DM from Nurse Interview.")
		self.dm['gdm_ni'] = pkall(pd.concat([pheno.findinfield('20002','1221',instances), pheno.sex == 'female'], axis='columns'), axis='columns')
		logger.info(f"Init: {self.dm['gdm_ni'].sum()} Subjects with Gestational DM from Nurse Interview.")
		self.dm['t1dm_ni'] = pheno.findinfield('20002','1222',instances)
		logger.info(f"Init: {self.dm['t1dm_ni'].sum()} Subjects with Type 1 DM from Nurse Interview.")
		self.dm['t2dm_ni'] = pheno.findinfield('20002','1223',instances)
		logger.info(f"Init: {self.dm['t2dm_ni'].sum()} Subjects with Type 2 DM from Nurse Interview.")
		self.dm['anynsgt1t2_ni'] = pd.Series(self.dm[['alldm_ni', 'gdm_ni', 't1dm_ni', 't2dm_ni']].any(axis='columns'), dtype='boolean')
		logger.info(f"Init: {self.dm['t2dm_ni'].sum()} Subjects with any DM from Nurse Interview.")

		# Fill dm with values - Medication
		self.dm['drug_ins_sr'] = pheno.findinfield(['6153','6177'],'3',instances)
		logger.info(f"Init: {self.dm['drug_ins_sr'].sum()} Subjects with Insulin, Medication from Touchscreen.")
		self.dm['insat1yr'] = pheno.findinfield('2986','1',instances)
		logger.info(f"Init: {self.dm['insat1yr'].sum()} subjects with Insulin started within 1 yr of diagnosis from Touchscreen.")
		self.dm['drug_ins_ni'] = pheno.findinfield('20003','1140883066',instances)
		logger.info(f"Init: {self.dm['drug_ins_ni'].sum()} Subjects with Insulin Product, Medication from Nurse Interview.")
		self.dm['drug_metf_ni'] = pheno.findinfield('20003',['1140884600','1140874686','1141189090'],instances)
		logger.info(f"Init: {self.dm['drug_metf_ni'].sum()} Subjects with Metformin, Medication from Nurse Interview.")

		Glitazones    = ['1141171646','1141171652','1141153254','1141177600','1141177606']
		Meglitinides  = ['1141173882','1141173786','1141168660']
		Sulfonylureas = ['1140874718','1140874744','1140874746','1141152590','1141156984','1140874646','1141157284','1140874652','1140874674','1140874728']
		OtherOAD      = ['1140868902','1140868908','1140857508']
		self.dm['drug_nonmetf_oad_ni'] = pheno.findinfield('20003', Glitazones + Meglitinides + Sulfonylureas + OtherOAD, instances)
		logger.info(f"Init: {self.dm['drug_nonmetf_oad_ni'].sum()} Subjects with Non-metformin oral anti-diabetic drug, Medication from Nurse Interview.")

		# Age at Diagnosis combined from TS and NI (Remember: .loc[] enforces index/column names so this works as intended)
		self.dm['agedm_ts_or_ni'] = pheno.findfield('2976_0').mean(axis='columns')                        # Touchscreen - gestational DM
		self.dm['agediag_gdm_ni'] = pheno.findinterpolated('20009','20002','1220').mean(axis='columns') # Nurse interview - gestational DM
		self.dm.loc[self.dm['alldm_ni'], 'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1220').mean(axis='columns') # Nurse interview - all DM
		self.dm.loc[self.dm['gdm_ni'],   'agedm_ts_or_ni'] = self.dm.loc[self.dm['gdm_ni'], 'agediag_gdm_ni']                    # Nurse interview - gestational DM
		self.dm.loc[self.dm['t1dm_ni'],  'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1222').mean(axis='columns') # Nurse interview - type 1 DM
		self.dm.loc[self.dm['t2dm_ni'],  'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1223').mean(axis='columns') # Nurse interview - type 2 DM
		logger.info(f"Init: {sum(self.dm['agedm_ts_or_ni'] > 0)} subjects with age at diagnosis.")
#		logger.debug(f"Init: Subject={self.DEBUGsubject} {self.dm.loc[self.DEBUGsubject,]}")

		# Precalculate the prevalence
		self.prevalence = self.prevalenceA()
		self.prevalence = self[self._prevalence == self.T1Moderate].prevalenceB()
		self.prevalence = self[self._prevalence == self.T2Moderate].prevalenceC()

	@Eastwood.baseline.setter
	def baseline(self, value):
		"""Overrides setter to do some asserts."""
		assert value > self.UKBstartdate, f"Prevalence calculations on UKBiobank data are too unreliable for baselines prior to the UKBiobank start date ({self.UKBstartdate.date()})"
		if value > self.UKBbaseline:
			logger.warning(f"Prevalence calculations on UKBiobank data is somewhat unreliable for baselines prior to the UKBiobank end of assessment date ({self.UKBbaseline.date()}).")
		self._baseline = value

	@property
	def prevalence(self):
		"""Implements prevalence algorithm from Eastwood2016 paper (See Fig 2 in paper)."""
		out = self._prevalence
		if self.styletr is not None:
			logger.info(f"Converting output to style: {self.style}")
			out = pd.Series(out.apply(self.styletr), dtype='category')
		return out

	@prevalence.setter
	def prevalence(self, value):
		"""Setter used to calculate the prevalence."""
		if hasattr(value, 'index'):
			self._prevalence[value.index] = value
		else:
			logger.error(f"Oops! You set prevalence to something fishy...")

	@property
	def style(self):
		"""Return: self._style. Needed for setter."""
		return self._style

# Don't much like setting the translator here, but it's meaningless to change the style without changing the translator
#	The clean way would be to make a style class. A thought for future update?
	@style.setter
	def style(self, value):
		"""Assert value and set it to self._style if valid."""
		value = value.lower()
		assert value in [s.lower() for s in self.styles], f"Given style {value} is unknown. Supported styles include {self.styles}"
		self._style = value
		if self.style == 'eastwood':
			self.styletr = None
		elif self.style == 't1d':
			self.styletr = lambda x: self.binaryPositive if x in [self.T1Moderate, self.T1High] else (self.binaryNegative if x in [self.Negative] else self.binaryOther)
		elif self.style == 't2d':
			self.styletr = lambda x: self.binaryPositive if x in [self.T2Moderate, self.T2High] else (self.binaryNegative if x in [self.Negative] else self.binaryOther)

	def prevalenceA(self):
		"""Prevalence algorithm A - Distinguishes between diabetes presence/absence 

		# Prevalence 1.1:
		Identifies participants in whom diabetes unlikely.

		Based on absence of self-reported diabetes diagnoses or medication. The exception was diabetes
		(non-GDM) report by TS alone (n=1,025), which was not used to indicate possible diabetes, since only
		3% (36/1,025) reported other markers of diabetes, compared with 73% (18,572/25,383) of those with
		diabetes diagnoses by TS and NI, and 16% (35/217) of those with NI diagnoses alone.

		Presence of diabetes complications was not used to rule in/ out cases, due to low prevalence (0.3%).

		# Prevalence 1.2:
		Identifies possible GDM only.

		Based on TS or NI self-report of GDM. Participants excluded at this step did not have any evidence of
		current diabetes (medication or self-report of other types of diabetes), i.e. those who may have had
		pre-existing diabetes or developed type 1 or 2 diabetes subsequent to GDM continue in the algorithm.
		This included a total of 338/1072 (32%) women reporting GDM by TS and 41/234 (18%) by NI with evidence
		of current diabetes who were not currently pregnant. In addition, 22 women who reported age at GDM
		diagnosis of >=50 years also continued in the algorithm (i.e. were not assigned possible GDM status).

		# Prevalence 1.3:
		Identifies possible type 2 diabetes from use of non-metformin oral anti-diabetic drugs.

		Use of this type of medication likely to be exclusive to type 2 diabetes.

		# Prevalence 1.4:
		Identifies possible type 2 diabetes from ethnicity-specific age at diagnosis.

		Age of onset of diabetes is younger in the UK's largest ethnic minority groups; South Asians and
		African Caribbeans. We applied a priori-determined ethnic-specific cut-points to define older vs.
		younger age at diagnosis - indicating higher likelihood of type 2 vs. type 1 diagnoses respectively.
		Age of onset was taken from NI data if NI and TS data were available, and taken from TS data if NI data
		were not available, 1243 participants with no age data went into the "possible type 2 diabetes" arm.

		# Prevalence 1.5:
		Identifies possible type 1 diabetes from either insulin use (TS or NI), insulin commencement within
		1 year of diagnosis, or self-report of type 1 diabetes (NI).

		Almost all participants with type 1 diabetes will be on insulin; coupled with a lack of non-metformin
		oral anti-diabetic drugs and a younger age at onset (see rules 1.4 and 1.5), specificity for type 1
		diabetes is increased.

		Immediate commencement of insulin after diagnosis is obligatory for type 1diabetes. We included the
		criterion of commencement of insulin at less than 12 months post-diagnosis in addition to current
		insulin use as 96 people entering rule 1.6 reported insulin use at less than 12 months but no current
		insulin. These people are then carried forward to the next flowchart for finalising type 1 diabetes
		status, and are labelled as "possible type 1 diabetes". Otherwise they would be excluded, but it is
		possible to have been on insulin initially and not currently, e.g. due to pancreatic transplant/
		incorrect diagnosis.

		Characteristics of participants self-reporting (NI) type 1 diabetes at this stage in the algorithm
		appear to be congruent with the diagnosis +/- their average age at diagnosis was 21±9 years, 98% (TS
		273/278) or 95% (NI, 265/278) were on current insulin, only 10% (30/278) were on metformin and 92%
		(194/211) started insulin within 1 year of diagnosis. However, due to the low number of participants
		specifying this diagnosis at NI overall 0.09% (428/502,665), we sought other evidence of type 1
		diabetes in addition. NB: there was no option to specify diabetes type on TS self-report.
		"""

		prevalence = self._prevalence[self.dm.index]
		subjects_left = pd.Series(True, index=self.dm.index, dtype='boolean')
		logger.info(f"Prevalence Algorithm A Starting: Analysing {subjects_left.sum()} subjects.")

		# Fig 2 (Flowchart): 1.1
		x = self.dm[['gdmonly_sr', 'alldm_ni', 'gdm_ni', 't1dm_ni', 't2dm_ni', 'drug_ins_ni', 'drug_ins_sr', 'drug_metf_ni', 'drug_nonmetf_oad_ni']].any(axis='columns')
		# Next line diverges from Eastwood who does not consider in-patient data prior to the UKB baseline
		x = pd.concat([x, self.dm['date_anydm_ip'] <= self.baseline], axis='columns').any(axis='columns')
		prevalence[x != True] = self.Negative
		subjects_left[x != True] = False
		logger.info(f"   Prevalence 1.1: {sum(x != True)} subjects assigned '{self.Negative}'; {subjects_left.sum()} subjects remaining.")
#		logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")

		# Fig 2 (Flowchart): 1.2
		self.dm['anydmrx_ni_sr'] = self.dm[['drug_ins_ni', 'drug_metf_ni', 'drug_nonmetf_oad_ni', 'drug_ins_sr']].any(axis='columns')

		x = pd.DataFrame({1: pkall(pd.DataFrame({1: self.dm['gdmonly_sr'],
                                                         2: self.dm[['anydmrx_ni_sr', 't1dm_ni', 't2dm_ni']].any(axis='columns') != True}),
		                           axis='columns'),
		                  2: pkall(pd.DataFrame({1: self.dm['gdm_ni'],
                                                         2: self.dm['agediag_gdm_ni'] < 50,
		                                         3: self.dm[['anydmrx_ni_sr', 't1dm_ni', 't2dm_ni']].any(axis='columns') != True}),
		                           axis='columns')
		    }).any(axis='columns')
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.GDModerate
		subjects_left[x] = False
		logger.info(f"   Prevalence 1.2: {x.sum()} subjects asigned '{self.GDModerate}'; {subjects_left.sum()} subjects remaining.")
#		logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")

		# Fig 2 (Flowchart): 1.3
		x = self.dm['drug_nonmetf_oad_ni']
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.T2Moderate
		subjects_left[x] = False
		logger.info(f"   Prevalence 1.3: {x.sum()} subjects asigned '{self.T2Moderate}'; {subjects_left.sum()} subjects remaining.")
#		logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")

		# Fig 2 (Flowchart): 1.4
		x = pd.concat([
		        self.dm['agedm_ts_or_ni'] > 36,
		        pkall(pd.concat([self.dm['agedm_ts_or_ni'] >= 31 , self.dm['ethnic_sa_afc']], axis='columns'), axis='columns'),
		    ], axis='columns').any(axis='columns')
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.T2Moderate
		subjects_left[x] = False
		logger.info(f"   Prevalence 1.4: {x.sum()} subjects asigned '{self.T2Moderate }'; {subjects_left.sum()} subjects remaining.")
#		logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")

		# Fig 2 (Flowchart): 1.5
		x = self.dm[['drug_ins_sr', 'drug_ins_ni', 'insat1yr', 't1dm_ni']].any(axis='columns')
		# Next line diverges from Eastwood who does not consider in-patient data prior to the UKB baseline
		x = pd.concat([x, self.dm['date_t1dm_ip'] <= self.baseline], axis='columns').any(axis='columns')
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		subjects_left[x] = False
		prevalence.loc[subjects_left != False] = self.T2Moderate
		logger.info(f"   Prevalence 1.5: {sum(subjects_left != False)} subjects asigned '{self.T2Moderate}'; {x.sum()} subjects remaining.")
#		logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		prevalence[x] = self.T1Moderate
		logger.info(f"Prevalence Algorithm A finished: Remaining {x.sum()} subjects asigned '{self.T1Moderate}'.")

		return prevalence[self.dm.index]

	def prevalenceB(self):
		"""Prevalence algorithm B. Classify type 1 diabetes into probable or possible.

		Prevalence 2.1:
		Identifies probable type 1 diabetes from self-report of type 1 diabetes, nurse interview

		Deemed to indicate probable type 1 diabetes - see above. NB - 10 participants co-reported type 1 and
		type 2 diabetes at NI, these were re-classified as non-specific diabetes.

		Prevalence 2.2:
		Identifies probable type 1 diabetes from insulin commencement within 1 year of diagnosis and current
		insulin use.

		Of those subjects entering this algorithm without nurse report of type 1 diabetes, some received
		insulin within a year of diagnosis, some self reported current insulin use on touchscreen, and some
		reported current insulin use to the nurse.  Metformin use was reported in some of those with "possible"
		versus those with "probable" type 1 diabetes.
		"""
		prevalence = self._prevalence[self.dm.index]
		subjects_left = pd.Series(True, index=self.dm.index, dtype='boolean')
		logger.info(f"Prevalence Algorithm B Starting: Analysing {subjects_left.sum()} subjects.")

		# Fig 2 (Flowchart): 2.1
		x = self.dm['t1dm_ni']
		# Next line diverges from Eastwood who does not consider in-patient data prior to the UKB baseline
		x = pd.concat([x, self.dm['date_t1dm_ip'] <= self.baseline], axis='columns').any(axis='columns')
		prevalence[x] = self.T1High
		subjects_left[x] = False
		logger.info(f"   Prevalence 2.1: {x.sum()} subjects asigned '{self.T1High}'; {sum(subjects_left != False)} subjects remaining.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm B.")

		# Fig 2 (Flowchart): 2.2
		x = pkall(pd.concat([self.dm['insat1yr'], self.dm[['drug_ins_sr', 'drug_ins_ni']].any(axis='columns')], axis='columns'), axis='columns')
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.T1High
		subjects_left[x] = False

		logger.info(f"   Prevalence 2.2: {x.sum()} subjects asigned '{self.T1High}'; {subjects_left.sum()} subjects remaining.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm B.")

		prevalence[subjects_left] = self.T1Moderate
		logger.info(f"Prevalence Algorithm B finished: Remaining {subjects_left.sum()} subjects asigned '{self.T1Moderate}'.")
		return prevalence[self.dm.index]

	def prevalenceC(self):
		"""Prevalence algorithm C. Classify type 2 diabetes into probable or possible.

		3.1
		Identifies participants whose only diabetic medication is metformin

		Metformin can be used to treat conditions other than diabetes, e.g. polycystic ovarian syndrome,
		therefore treatment with metformin only is less specific for diabetes than if used in conjunction with
		other anti-diabetic drugs.

		3.2
		Identifies participants whose only diabetic medication is metformin but did not report diabetes at
		nurse interview

		With the exception of 12/169 (7%) who reported diabetes on TS, no other diabetes self-report data was
		present for participants who were deemed unlikely to have diabetes at this stage.

		3.3
		Identifies participants who are on oral anti-diabetic agents other than metformin (i.e. likely to be
		specific for type 2 diabetes)

		Participants assigned "probable type 2 diabetes" status at this stage were more likely than those who
		continued in the algorithm to have developed diabetes at an older age (55±9 vs. 52±9 years, p<0.00
		and to be taking metformin concurrently (79%, 5404/6809 vs. 51%, 8750/17152, p<0.001), and less likely
		to be on insulin (8%, 546/6809 vs. 51%, 8750/17152, p<0.001).

		3.4
		Identifies participants not on insulin (i.e. unlikely to have type 1 diabetes)

		Participants who were assigned "probable type 2 diabetes" status (i.e. were not receiving insulin) at
		this stage were more likely than those who continued in the algorithm to have developed diabetes at an
		older age (56±9 vs. 46±9 years, p<0.001), though rates of metformin use were simil (51%, 6983/13761
		vs.  52%, 1767/3391, p=0.2).

		3.5
		Identifies probable type 1 diabetes on the basis of self-report of type 1 diabetes, nurse interview

		Participants who were assigned "probable type 1 diabetes" status at this stage were more likely than
		those who were assigned "possible type 2 diabetes" status to have developed diabetes at a younger age
		(48±8 vs. 50±7 years, p<0.1) and to have commenced insulin at less than one year after diagnis
		(72%, n=79/109 vs. 39%, n=1197/3040, p<0.001), and less likely to be taking metformin currently (27%,
		n=32/118 vs. 53%, n=1735/3273, p<0.001).

		Participants were assigned "possible type 2 diabetes" rather than "probable type 2 diabetes" status at
		this stage as insulin use and relatively high rates of commencement of insulin within a year of
		diagnosis (39%, n=1197/3040) cast doubt over a type 2 diagnosis.
		"""
		prevalence = self._prevalence[self.dm.index]
		subjects_left = pd.Series(True, index=self.dm.index, dtype='boolean')
		logger.info(f"Prevalence Algorithm C Starting: Analysing {subjects_left.sum()} subjects.")

		# Fig 2 (Flowchart): 3.1
		x = pkall(pd.concat([self.dm['drug_metf_ni'],
		                     self.dm[['drug_ins_sr', 'drug_ins_ni', 'drug_nonmetf_oad_ni']].any(axis='columns') != True], axis='columns'),
		          axis='columns')
		logger.info(f"   Prevalence 3.1: {x.sum()} subjects directed to 3.2; {x.size - x.sum()} subjects directed to 3.3.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm C.")

		# Fig 2 (Flowchart): 3.2
		y = pkall(pd.concat([x, self.dm['anynsgt1t2_ni'] != True], axis='columns'), axis='columns')
		# Next line diverges from Eastwood who does not consider in-patient data prior to the UKB baseline
		y[self.dm['date_anydm_ip'] <= self.baseline] = False
		prevalence[y] = self.Negative
		logger.info(f"   Prevalence 3.2: {y.sum()} subjects asigned '{self.Negative}'; remaining {x.sum() - y.sum()} subjects directed to 3.3.")
		subjects_left[y] = False
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm C.")

		# Fig 2 (Flowchart): 3.3
		x = self.dm['drug_nonmetf_oad_ni']
		# Next line diverges from Eastwood who does not consider in-patient data prior to the UKB baseline
		x = pd.concat([x, self.dm['date_t2dm_ip'] <= self.baseline], axis='columns').any(axis='columns')
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.T2High
		subjects_left[x] = False
		logger.info(f"   Prevalence 3.3: {x.sum()} subjects asigned '{self.T2High}'; {subjects_left.sum()} subjects remaining.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm C.")

		# Fig 2 (Flowchart): 3.4
		x = self.dm[['drug_ins_sr', 'drug_ins_ni']].any(axis='columns') != True
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		prevalence[x] = self.T2High
		subjects_left[x] = False
		logger.info(f"   Prevalence 3.4: {x.sum()} subjects asigned '{self.T2High}'; {subjects_left.sum()} subjects remaining.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm C.")

		# Fig 2 (Flowchart): 3.5
		x = self.dm['t1dm_ni']
		x = pkall(pd.concat([x, subjects_left], axis='columns'), axis='columns')
		subjects_left[x] = False
		prevalence[subjects_left] = self.T2Moderate

		logger.info(f"   Prevalence 3.5: {subjects_left.sum()} subjects asigned '{self.T2Moderate}'; {x.sum()} subjects remaining.")
		try: logger.debug(f"     Testsubject: {self.DEBUGsubject} = '{prevalence[self.DEBUGsubject]}'")
		except KeyError:
			logger.debug(f"     Testsubject: '{self.DEBUGsubject}' not in Algorithm C.")

		prevalence[x] = self.T1High
		logger.info(f"Prevalence Algorithm C Finished: Remaining {x.sum()} subjects asigned '{self.T1High}'.")
		return prevalence[self.dm.index]

	def to_incidence(self, *args, **kwargs):
		"""Convert to Incidence class (with equal args, like baseline)."""
		out = Incidence(*args, prev=self, **kwargs)
		return out



###########################################################
#
# --%%  DEFINE: Incidence class from the Eastwood2016 paper.  %%--

class Incidence(Prevalence):
	"""Implements incidence algorithm 'B' from Eastwood2016 paper (See Fig 3 in paper).
	Only implements algorithm B since A requires CPRD data not available to us."""

	# UKB End date
	UKBenddate = pd.Timestamp.max

	def __init__(self, *args, enddate=None, interval=None, prev=None, **kwargs):
		"""Init the Incidence object."""
		if prev is None:
			super().__init__(*args, **kwargs)
		else:
			self = copy.deepcopy(prev)
		self.enddate = enddate
		self.interval = interval
		self.dm['prevalent'] = self._prevalence.isin([self.T1High, self.T2High])
		logger.info(f"Incicence.init: {self.dm['prevalent'].sum()} subjects with DM diagnosis prior to baseline ({self.baseline.date()})")

	@Eastwood.baseline.setter
	def baseline(self, value):
		"""Overrides setter to do some asserts."""
		assert value >= self.UKBbaseline, f"Incidence calculations on UKBiobank data are too unreliable for baselines prior to the UKBiobank end of assessment date ({self.UKBbaseline.date()})."
		self._baseline = value

	@property
	def enddate(self):
		"""Getter for enddate. Needed for the setter."""
		return self._enddate

	@enddate.setter
	def enddate(self, value):
		"""Enddate Setter. Ensure that enddate is a datetime object."""
		if value is None:
			self._enddate = self.UKBenddate
		self._enddate = value if isinstance(value, datetime) else pd.to_datetime(value)

	@property
	def interval(self):
		"""Getter for interval. Needed for the setter."""
		return self._interval

	@interval.setter
	def interval(self, value):
		"""Interval Setter: Ensure that interval is a timedelta object."""
		if value is None:
			self._interval = None
		else:
			self._interval = pd.Timedelta(value)

	@property
	def anydm(self):
		"""Incidence of Type-1 + Type-2 + Unspecified Diabetes Mellitus."""
		incidence = self._incidence[self.dm.index]
		x = pkall(pd.concat([~self.dm['prevalent'], self.dm['date_anydm_ip'] > self.baseline, self.dm['date_anydm_ip'] < self.enddate], axis='columns'), axis='columns')
		incidence[x] = self.dm.loc[x, 'date_anydm_ip']
		logger.info(f"Incidence Any DM: {x.sum()} subjects with Any Diabets diagnosis data.")
		if self.interval:
			incidence = (incidence - self.baseline).floordiv(self.interval) + 1
		return incidence[self.dm.index]

	@property
	def t1dm(self):
		"""Incidence of Type-1 Diabetes Mellitus."""
		incidence = self._incidence[self.dm.index]
		x = pkall(pd.concat([~self.dm['prevalent'], self.dm['date_t1dm_ip'] > self.baseline, self.dm['date_t1dm_ip'] < self.enddate], axis='columns'), axis='columns')
		incidence[x] = self.dm.loc[x, 'date_t1dm_ip']
		logger.info(f"Incidence Type-1 DM: {x.sum()} subjects with Type-1 DM diagnosis data.")
		if self.interval:
			incidence = (incidence - self.baseline).floordiv(self.interval) + 1
		return incidence[self.dm.index]

	@property
	def t2dm(self):
		"""Incidence of Type-2 Diabetes Mellitus."""
		incidence = self._incidence[self.dm.index]
		x = pkall(pd.concat([~self.dm['prevalent'], self.dm['date_t2dm_ip'] > self.baseline, self.dm['date_t2dm_ip'] < self.enddate], axis='columns'), axis='columns')
		incidence[x] = self.dm.loc[x, 'date_t2dm_ip']
		logger.info(f"Incidence Type-2 DM: {x.sum()} subjects with Type-2 DM diagnosis data.")
		if self.interval:
			incidence = (incidence - self.baseline).floordiv(self.interval) + 1
		return incidence[self.dm.index]

	def to_incidence(self, *args, **kwargs):
		"""Dummy converter."""
		return self

