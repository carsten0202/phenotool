

###########################################################
#
# --%%  Setup and Initialize  %%--

import copy
import logging
import pandas as pd
import sys

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)

#import phenotool.options as OPTIONS
#from pklib.pkclick import CSV, gzFile




###########################################################
#
# --%%  DEFINE: Eastowwd class from the Eastwood2016 paper.  %%--

class Eastwood():
	"""Implements a Diabetes classification based on Eastwood2016."""

	# Define our categories
	Negative   = "Diabetes_Unlikely"
	T1High     = "Probable_Type1"
	T2High     = "Probable_Type2"
	T1Moderate = "Possible_Type1"
	T2Moderate = "Possible_Type2"
	GDModerate = "Possible_Gestational"

	def __init__(self, pheno, agediag=None, ethnicity=None, high=None, moderate=None, treatments=None):
		"""Based on the Eastwood2016 paper.

		The following refers to the prevalence algorithm only.
		"""

		# Start with some assessments. Do we have the data?
		#	Remember that we mau want to use this definition on other datasets; not just UKB.
		#	This is actually an argument for making Eastwood a command and not an option...
		#		It should then be above the output format commands; eg: phenotool eastwood snptest
		#		Will make assessment critical, since data could be missing (and user will need to provide much information...)

		# Possible output (Categorical):
			# A missing value should indicate 'unassessed' since that's an option if input data are missing.

		# Initialize attributes
		self.dm = pd.DataFrame(index=pheno.index)
		self.subject_index = pheno._obj.index 
		self._prevalence = pd.Series(
		    pd.Categorical([pd.NA] * pheno.index.size,
		                   categories=[self.Negative, self.T1High, self.T2High, self.T1Moderate, self.T2Moderate, self.GDModerate]),
		    name="Prevalence",
		    index=pheno.index
		) # This guy is a frozen init; actual prevalence is returned through a property getter


		# Fill dm with values - Touchscreen
		self.dm['gdmonly_sr'] = pd.concat([pheno.findinfield('4041','1'), pheno.sex == 'female'], axis='columns').all(axis='columns')
		logger.info(f"Eastwood: {self.dm['gdmonly_sr'].sum()} subjects with Gestational Diabetes from Touchscreen.")

		# Fill dm with values - Nurse Interview
		self.dm['alldm_ni'] = pheno.findinfield('20002','1220')
		logger.info(f"Eastwood: {self.dm['alldm_ni'].sum()} Subjects with any type DM from Nurse Interview.")
		self.dm['gdm_ni'] = pd.concat([pheno.findinfield('20002','1221'), pheno.sex == 'female'], axis='columns').all(axis='columns')
		logger.info(f"Eastwood: {self.dm['gdm_ni'].sum()} Subjects with Gestational DM from Nurse Interview.")
		self.dm['t1dm_ni'] = pheno.findinfield('20002','1222')
		logger.info(f"Eastwood: {self.dm['t1dm_ni'].sum()} Subjects with Type 1 DM from Nurse Interview.")
		self.dm['t2dm_ni'] = pheno.findinfield('20002','1223')
		logger.info(f"Eastwood: {self.dm['t2dm_ni'].sum()} Subjects with Type 2 DM from Nurse Interview.")

		# Fill dm with values - Medication
		self.dm['drug_ins_sr'] = pheno.findinfield(['6153','6177'],'3')
		logger.info(f"Eastwood: {self.dm['drug_ins_sr'].sum()} Subjects with Insulin, Medication from Touchscreen.")
		self.dm['insat1yr'] = pheno.findinfield('2986', '1')
		logger.info(f"Eastwood: {self.dm['insat1yr'].sum()} subjects with Insulin started within 1 yr of diagnosis from Touchscreen.")
		self.dm['drug_ins_ni'] = pheno.findinfield('20003','1140883066')
		logger.info(f"Eastwood: {self.dm['drug_ins_ni'].sum()} Subjects with Insulin Product, Medication from Nurse Interview.")
		self.dm['drug_metf_ni'] = pheno.findinfield('20003',['1140884600','1140874686','1141189090'])
		logger.info(f"Eastwood: {self.dm['drug_metf_ni'].sum()} Subjects with Metformin, Medication from Nurse Interview.")

		Glitazones    = ['1141171646','1141171652','1141153254','1141177600','1141177606']
		Meglitinides  = ['1141173882','1141173786','1141168660']
		Sulfonylureas = ['1140874718','1140874744','1140874746','1141152590','1141156984','1140874646','1141157284','1140874652','1140874674','1140874728']
		OtherOAD      = ['1140868902','1140868908','1140857508']
		self.dm['drug_nonmetf_oad_ni'] = pheno.findinfield('20003', Glitazones + Meglitinides + Sulfonylureas + OtherOAD)
		logger.info(f"Eastwood: {self.dm['drug_nonmetf_oad_ni'].sum()} Subjects with Non-metformin oral anti-diabetic drug, Medication from Nurse Interview.")

		# Age at Diagnosis combined from TS and NI (Remember: .loc[] enforces index/column names so this works as intended)
		self.dm['agedm_ts_or_ni'] = pheno.findfield('2976').mean(axis='columns')   # Touchscreen - gestational DM
		self.dm['agediag_gdm_ni'] = pheno.findinterpolated('20009','20002','1220') # Nurse interview - gestational DM
		self.dm.loc[self.dm['alldm_ni'], 'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1220')   # Nurse interview - all DM
		self.dm.loc[self.dm['gdm_ni'],   'agedm_ts_or_ni'] = self.dm.loc[self.dm['gdm_ni'], 'agediag_gdm_ni'] # Nurse interview - gestational DM
		self.dm.loc[self.dm['t1dm_ni'],  'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1222')   # Nurse interview - type 1 DM
		self.dm.loc[self.dm['t2dm_ni'],  'agedm_ts_or_ni'] = pheno.findinterpolated('20009','20002','1223')   # Nurse interview - type 2 DM
		logger.info(f"Eastwood: {sum(self.dm['agedm_ts_or_ni'] > 0)} subjects with age at diagnosis.")

		# Fill dm with values - Main ethnic groups (NB: this makes some assumptions i.e. British is white, and does not include mixed, Chinese or other Asian)
		self.dm.loc[pheno.findinfield('21000', ['1', '1001', '1002', '1003']), 'ethnic'] = 1 # White European
		self.dm.loc[pheno.findinfield('21000', ['3', '3001', '3002', '3003']), 'ethnic'] = 2 # South Asian
		self.dm.loc[pheno.findinfield('21000', ['4', '4001', '4002', '4003']), 'ethnic'] = 3 # African Caribbean
		self.dm.loc[pheno.findinfield('21000', ['2', '2001', '2002', '2003', '2004', '5', '6']), 'ethnic'] = 4 # Mixed or Other
		self.dm['ethnic_sa_afc'] = [ethnic in [2, 3] for ethnic in self.dm['ethnic']] # SA and AFC vs all other variable
		logger.info(f"Eastwood: Ethnic breakdown = {dict(zip(self.dm['ethnic'].value_counts().sort_index(), ['White European', 'South Asian', 'African Caribbean', 'Mixed or Other']))}.")

		# Validation
		Eastwood._validate(self)

	def __getitem__(self, key):
		"""Propagates index operations across the relevant attributes."""
		out = copy.deepcopy(self)
		out.dm = self.dm.loc[key,]
		out.subject_index = self.subject_index[key]
		out._prevalence = self._prevalence[key]
		return out

	@staticmethod
	def _validate(self):
		"""Check stuff. But what?

		If a parameter matches all or none of the subjects?
		"""
		pass

	@property
	def prevalence(self):
		"""Convienience wrapper which implements algorithms A+B+C"""
		out = self.prevalenceA()
		t1dm = out == self.T1Moderate
		t2dm = out == self.T2Moderate
		out[t1dm] = self[t1dm].prevalenceB()
		out[t2dm] = self[t2dm].prevalenceC()
		return out

	def prevalenceA(self):
		"""Implements prevalence algorithm A from Eastwood2016 paper (See Fig 2).

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
		subjects_left = pd.Series(True, index=self.dm.index) # Some tracker thingy; subjects_left? remaining?

		# Fig 2 (Flowchart): 1.1
		x = self.dm[['gdmonly_sr', 'alldm_ni', 'gdm_ni', 't1dm_ni', 't2dm_ni', 'drug_ins_ni', 'drug_ins_sr', 'drug_metf_ni', 'drug_nonmetf_oad_ni']].any(axis='columns')
		prevalence[~x] = self.Negative
		subjects_left[~x] = False
		logger.info(f"Prevalence 1.1: {sum(~x)} subjects assigned '{self.Negative}'; {sum(subjects_left != False)} subjects remaining.")
#		if not subjects_left.any():
#			sys.exit("No more subjects left. Need stop code here.")

		# Fig 2 (Flowchart): 1.2
		self.dm['anydmrx_ni_sr'] = self.dm[['drug_ins_ni', 'drug_metf_ni', 'drug_nonmetf_oad_ni', 'drug_ins_sr']].any(axis='columns')
		x = pd.DataFrame({1: pd.DataFrame({1: self.dm['gdmonly_sr'],
                                                   2: ~self.dm[['anydmrx_ni_sr', 't1dm_ni', 't2dm_ni']].any(axis='columns'),
		                     }).all(axis='columns'),
		                  2: pd.DataFrame({1: self.dm['gdm_ni'],
                                                   2: self.dm['agediag_gdm_ni'] < 50,
		                                   3: ~self.dm[['anydmrx_ni_sr', 't1dm_ni', 't2dm_ni']].any(axis='columns'),
		                     }).all(axis='columns')
		    }).any(axis='columns')
		x = pd.concat([x, subjects_left], axis='columns').all(axis='columns')
		prevalence[x] = self.GDModerate
		subjects_left[x] = False
		logger.info(f"Prevalence 1.2: {x.sum()} subjects asigned '{self.GDModerate}'; {sum(subjects_left != False)} subjects remaining.")

		# Fig 2 (Flowchart): 1.3
		x = self.dm['drug_nonmetf_oad_ni']
		x = pd.concat([x, subjects_left], axis='columns').all(axis='columns')
		prevalence[x] = self.T2Moderate
		subjects_left[x] = False
		logger.info(f"Prevalence 1.3: {x.sum()} subjects asigned '{self.T2Moderate}'; {sum(subjects_left != False)} subjects remaining.")

		# Fig 2 (Flowchart): 1.4
		x = pd.concat([self.dm['agedm_ts_or_ni'] > 36,
		              pd.concat([self.dm['agedm_ts_or_ni'] >= 31 , self.dm['ethnic_sa_afc']], axis='columns').all(axis='columns'),
		    ], axis='columns').any(axis='columns')
		x = pd.concat([x, subjects_left], axis='columns').all(axis='columns')
		prevalence[x] = self.T2Moderate
		subjects_left[x] = False
		logger.info(f"Prevalence 1.4: {x.sum()} subjects asigned '{self.T2Moderate }'; {sum(subjects_left != False)} subjects remaining.")

		# Fig 2 (Flowchart): 1.5
		x = self.dm[['drug_ins_sr', 'drug_ins_ni', 'insat1yr', 't1dm_ni']].any(axis='columns')
		x = pd.concat([x, subjects_left], axis='columns').all(axis='columns')
		subjects_left[x] = False
		prevalence.loc[subjects_left != False] = self.T2Moderate
		logger.info(f"Prevalence 1.5: {sum(subjects_left != False)} subjects asigned '{self.T2Moderate}'; {x.sum()} subjects remaining.")
		prevalence[x] = self.T1Moderate
		logger.info(f"Prevalence Algorithm A finished: remaining {x.sum()} subjects asigned '{self.T1Moderate}'.")

		# Report and return
		return prevalence[self.subject_index]

	def prevalenceB(self):
		"""Implements prevalence algorithm B from Eastwood2016 paper (See Fig 2)."""
		prevalence = self._prevalence[self.dm.index]
		subjects_left = pd.Series(True, index=self.dm.index) # Some tracker thingy; subjects_left? remaining?
		# Fig 2 (Flowchart): 2.1
		# Fig 2 (Flowchart): 2.2

		prevalence[self.dm.index] = self.T1Moderate # Reset for testing
		return prevalence[self.subject_index]

	def prevalenceC(self):
		"""Implements prevalence algorithm C from Eastwood2016 paper (See Fig 2)."""
		prevalence = self._prevalence[self.dm.index]
		subjects_left = pd.Series(True, index=self.dm.index) # Some tracker thingy; subjects_left? remaining?
		# Fig 2 (Flowchart): 3.1
		# Fig 2 (Flowchart): 3.2
		# Fig 2 (Flowchart): 3.3
		# Fig 2 (Flowchart): 3.4
		# Fig 2 (Flowchart): 3.5

		prevalence[self.dm.index] = self.T2Moderate # Reset for testing
		return prevalence[self.subject_index]


