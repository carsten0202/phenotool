
#
# --%% plclick.py  %%--
#
 
import click
import logging




# TODO: Here's some fun. All click type extenders must obey thede five guys:
#    it needs a name
#    it needs to pass through None unchanged
#    it needs to convert from a string
#    it needs to convert its result type through unchanged (eg: needs to be idempotent)
#    it needs to be able to deal with param and context being None. This can be the case when the object is used with prompt inputs.


#
# -%  Class pkclick.gzFile  %-

class gzFile(click.File):
	"""A Class for detecting compressed files and automagically decrompress them."""
	magic_dict = {
		b"\x1f\x8b\x08"     : "gz",
		b"\x42\x5a\x68"     : "bz2",
		b"\x50\x4b\x03\x04" : "zip"
	}

	def convert(self, value, param, ctx):
		"""Converts (extracts) compressed input."""
		import gzip
		import io
		f = super().convert(value, param, ctx)
		try:
			if self._getziptype(f, self.magic_dict) is not None:
				return io.TextIOWrapper(gzip.GzipFile(fileobj=f))
		except UnicodeDecodeError:
			self.fail("Could not interpret input. Did you remember to use binary mode? eg gzFile(mode='rb')")
		return io.TextIOWrapper(f)

	@staticmethod
	def _getziptype(f, magic_dict):
		"""If f is seekable, return the zip flavor."""
		if f.seekable():
			file_start = f.read(max(len(x) for x in magic_dict))
			f.seek(0)
			for magic, filetype in magic_dict.items():
				if file_start.startswith(magic):
					return filetype
		return None



#
# -%  CLASS: pkclick.csv  %-

class CSV(click.ParamType):
	"""A Class for handling option values that are comma separated values."""
	name = "csv_list"
	def convert(self, value, param, ctx):
		import csv
		value = super().convert(value, param, ctx)
		return next(csv.reader([value]))



#
# -%  CLASS: Samples  %-

class SampleList(gzFile):
	"""Obtain a list of samples from a file (or '-'... maybe?)."""
	def convert(self, value, param, ctx):
		f = super().convert(value, param, ctx)
		logging.debug(f"SampleList: Input is seekable: {f.seekable()}.")
		if self._isVCF(f):
			logging.info("SampleList: Treating samples file as VCF.")
			try:
				import vcf
				vcf_r = vcf.Reader(f, compressed=False)
				record = next(vcf_r)
				samples = [subject.sample for subject in record.samples]
			except ImportError:
				self.fail("Encountered VCF imput but could not find the PyVCF module. Either install with 'pip install PyVCF' or provide different input.")
		elif self._isTable(f):
			import pkcsv as csv
			(riter, hdr) = csv.reader(f)
			logging.info(f"SampleList: Treating samples file as Table/CSV. Reading from column '{hdr[0]}'.")
			samples = [row[0] for row in riter]
		else:
			logging.info("SampleList: Treating samples file as list of plain IDs.")
			samples = [line.rstrip() for line in f.readlines()]
		logging.debug(f"SampleList: samples[:5]={samples[:5]}.")
		logging.info(f"SampleList: Read {len(samples)} sample identifiers.")
		return samples
#		self.fail("Could not parse input. Please try a different file format.")

	@staticmethod
	def _isVCF(f):
		"""Is f a VCF file? Returns 'None' if it couldn't check f."""
		if f.seekable():
			file_start = f.readline()
			f.seek(0)
			return file_start.startswith("##fileformat=VCFv4")
		return None

	@staticmethod
	def _isTable(f):
		"""Is f a tsv/csv file? Returns 'None' if it couldn't check f."""
		if f.seekable():
			import csv
			file_start = f.readline()
			f.seek(0)
			dialect = csv.Sniffer().sniff(file_start)
			return len(next(csv.reader([file_start], dialect))) > 1
		return None

