###########################################################
#
# --%%  Outputting customizable text files with Phenotypes  %%--

import logging
import pandas as pd
import sys

from .pkpheno import Phenotype

assert sys.version_info >= (3, 8), f"{sys.argv[0]} requires Python 3.8.0 or newer. Your version appears to be: '{sys.version}'."
logger = logging.getLogger(__name__)



###########################################################
#
# --%%  DEFINE: TestFile class .  %%--

class TextFile(Phenotype):
	"""Holds phenotypes in a custom text format for flexible output.
	"""

	__name__ = "TextFile"

	MAGIC_COLS = {'ID'    : Phenotype.MAGIC_COLS['IID'],
	              'SEX'   : Phenotype.MAGIC_COLS['SEX']}

	mkey_id = "ID" # Also the index, so must be unique.

	def __init__(self, *args, **kwargs):
		"""Init the TextFile object."""
		super().__init__(*args, **kwargs)

