
#
# -%  Define shared help strings  %-

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

samples="""
File with samples to include in the output. Samples will be outputted in the exact same order as in the sample file
including outputting samples with missing values if no phenotype information was found. The sample file can be a plain text file with sample names or a VCF file with sample genotypes.
"""


