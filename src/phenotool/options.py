
#
# --%%  Define shared help strings  %%--

columns = """
Comma separated list of columns to output in addition to mandatory columns. Default is to output all columns.
"""

covariates = """
Comma separated list of columns with covariates. Print only these columns (plus any mandatory columns)
"""

fam = """
Output fam format for plink1.9. If specified, takes one mandatory argument stating the column name of the phenotype to
include in the sixth fam file column.
"""

files = """
Input File(s). %(prog)s accepts one or more input files including '-' symbolizing stdin. The precise format of each
input file will be autodetected, but it should be some form of delimited text data file like 'csv' or tab-delimited.
"""

log = """
Control logging. Valid levels: 'debug', 'info', 'warning', 'error', 'critical'.
"""

phenotypes = """
Comma separated list of columns with phenotypes.
"""

samples = """
File with samples to include in the output. Samples will be outputted in the exact same order as in the sample file
including outputting samples with missing values if no phenotype information was found. The sample file can be a plain
text file with sample names or a VCF file with sample genotypes.
"""



#
# -%  For Derive and friends  %-

columnname = """
Name of the output data column. Will be overridden if it already exists.
"""

columnprefix = """
Prefix for output data columns.
"""



#
# -%  For Snptest  %-

covariates = """
Comma separated list of columns with covariates. Print only these columns (plus any mandatory columns).
"""

phenotypes = """
Comma separated list of columns with phenotypes.
"""


#
# -%  For TextFile  %-

csv = """Sets output to comma-separated values.
"""

tsv = """Sets output to tab-separated values.
"""

#
# -%  For UKBioBank  %-

datafields = """
Data Field(s) to output. Several fields can be specified as a comma-separated string with no spaces. Mandatory.
"""

instances = """
Instances to output. Several instances can be specified as a comma-separated string with no spaces. Optional.
"""

values = """
Only subject having this value in at least one of the specified data field(s) will be printed. Note that this uses
text-based Regular Expression matching, so '-v 1' will match a value of '12'. Use '-v \"\\b1\\b\"' to match only '1'.
Default is to print all subjects having at least one non-empty value in the specified data field(s).
"""


#
# -%  For the Prevalence/Incidence Algorithms  %-

agediag = """
Subjects age at the time disease was reported.
"""

baseline = """
Baseline date. All prior information will be considered baseline data.
"""

datediag = """
Diagnoses and dates when they were first given.
"""

diagnosed = """
Empirical evidence eg. provided by Doctor.
"""

ethnicity = """
Used to set typical age of onset for type 2 diabetes.
"""

inciinterval = """
Give incidence as intervals of time. <TIME_STR> should be a timeunit or a number followed by timeunit (eg "day" or "5 days"). Suggested units: day, week, month and year.
"""

prevstyles = """
Sets the style of the prevalence output categories.
"""

reported = """
Circumstantial evidence eg. based on questionnaire, interview, etc.
"""

stopdate = """
Last date to consider. All information after this date will be ignored.
"""

treatments = """
Evidence of relevant treatments; eg. Insulin for Diabetes.
"""
