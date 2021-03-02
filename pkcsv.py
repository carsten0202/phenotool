
#
# --%%  pkcsv.py  %%--
#

# Import functions from csv that we want to make available upstairs

from csv import DictWriter, writer
import csv
import itertools
import logging

logger = logging.getLogger(__name__)

def reader(f, dialect=None, **fmtparams):
	"""Autodetects input with sniffer and returns tuple with iter and header."""
	(header, fout, dialect) = _configure(f)
	next(fout)
	riter = csv.reader(f, dialect=dialect, **fmtparams)
	return riter, header

def DictReader(f, fieldnames=None, dialect=None, *args, **kwds):
	"""Autodetects input with sniffer and returns dict with header as keys."""
	(header, fout, dialect) = _configure(f)
	next(fout)
	riter = csv.DictReader(fout, fieldnames=header, dialect=dialect, *args, **kwds)
	return riter

def _configure(f):
	"""Check input and find header and readable iterator."""
	(f1, f2) = itertools.tee(f,2)
	line1 = next(f1)
	line2 = next(f1)
	dialect = csv.Sniffer().sniff(line1 + line2)
	reader = csv.reader([line1, line2], dialect)
	header = next(reader)
	second = next(reader)
	logger.debug(f"{header} {second}")
	if len(header) + 1 == len(second):
		logger.info(f"R table format detected in input.")
		header = ["row.index"] + header
	return header, f2, dialect


