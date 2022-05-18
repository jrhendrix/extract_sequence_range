'''
FILE:	xRange.py
AUTHOR:	J.R. Hendrix
URL: 	http://stronglab.org
DESC:	This script extracts a range of nucleotides 
		from a fasta sequence file
'''

# IMPORT FROM PYTHON STANDARD LIBRARY
import argparse
import os
import sys

from Bio import SeqIO	# Source: https://biopython.org/wiki/SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


class File:
	""" Base class for all file-types """

	def __init__(self, path, file_type=None):
		self._path = None
		self.path = path
		self.file_type = file_type

	@property
	def path(self):
		return self._path
	
	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isfile(value):
			self._path = value
		else:
			raise FileNotFoundError(value)

	@property
	def filename(self):
		return os.path.basename(self.path)

	@property
	def file_prefix(self):
		return self.filename.split(".")[0]

	@property
	def extension(self):
		return self.filename.split(".")[-1]

class Dir:
	""" Base class for system directories """

	def __init__(self, path):
		self._path = None
		self.path = path

	@property
	def path(self):
		return self._path

	@path.setter
	def path(self, value):
		if not os.path.isabs(value):
			value = os.path.join(os.getcwd(), value)
		if os.path.isdir(value):
			self._path = value
		else:
			raise NotADirectoryError(value)

def modify_index(start, stop):

	# INDEX MODIFIER
	## Note: Input will be position. Translate to index
	start = start - 1
	stop = stop - 1

	# FIX OUT OF BOUNDS CASES
	if start < 0:
		start = 0
	if stop < 0:
		stop = 0

	# CHECK THAT 'STOP' COMES AFTER 'START'
	if stop < start: # flip indecies
		print('Start index is larger than stop index. Inverting...\n')
		iStart = start
		start = stop
		stop = iStart

	return(start, stop)




def extract_range(args, start, stop, outdir):
	'''Find a range of indexes and write bases to file'''

	start, stop = modify_index(start, stop)

	outname = ''.join((args.savename, '.fasta'))
	outfile = '/'.join((outdir.path, outname))
	f2 = open(outfile, 'w')

	# LOOP THROUGH SEQUENCE FILE
	records = []
	f1 = open(args.fasta, 'r')
	for record in SeqIO.parse(f1, 'fasta'):
		# Restrict extraction to a specific node
		if args.id != None and args.id != record.id:
			continue

		max_index = len(record.seq)

		# INDEX MODIFIER
		if stop > max_index:
			stop = max_index

		print('xRange uses an inclusive range')
		print(f'Exporting sequence at positions: {start+1} - {stop+1}')
		print(f'\t(i.e. index: {start} - {stop})\n')


		fragment = record.seq[start:stop+1]
		extract_length = len(fragment)
		print(f'Extract length: {extract_length}')

		xTract = '-'.join((str(start+1), str(stop+1)))
		xTract = ''.join(('extract:', xTract))

		if args.outID == None:
			newID = ' '.join((record.id, xTract))
		else:
			newID = ' '.join((args.outID, xTract))
		# TAKE REVERSE COMPLIMENT
		if args.reverse_complement:
			fragment = fragment.reverse_complement()
			newID = ' '.join((newID, 'reverse_complement'))


		new_record = SeqRecord(fragment, newID, '', '')
		records.append(new_record)

		break

	# WRITE RECORD(S) TO OUTPUT
	SeqIO.write(records, f2, "fasta")
	f1.close()
	f2.close()


def check_input(args):

	# CHECK IF FILE EXISTS
	try:
		f = File(args.fasta)
	except IOError:
		print("ERROR: Could not find file")
		return 1

	# CHECK THAT FILE IS FASTA FORMAT
	extensions = ('fasta', 'fa', 'fna', 'faa')
	if f.extension not in extensions:
		print("ERROR: Input file was not in FASTA format.")
		return 2

	return 0

def make_output(args):

	# CHECK THAT OUTPUT DIRECTORY EXISTS
	d = Dir(args.p)
	outpath = '/'.join((args.p, args.output_directory))
	if os.path.isdir(outpath):
		outdir = Dir(outpath)
	else:
		os.makedirs(outpath)
		outdir = Dir(outpath)

	return outdir



def main(program):
	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='xRange', add_help=False)
	parent_parser.add_argument('-f', '--fasta', help="FASTA file containing sequences")
	parent_parser.add_argument('--id', default=None, help='Match sequence ID')
	parent_parser.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parent_parser.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parent_parser.add_argument('-p', default=cwd, help='Path to output', type=str)
	parent_parser.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parent_parser.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	subparsers = parent_parser.add_subparsers(help='sub-command help')

	# PARSER : DEFINED
	def_parser = subparsers.add_parser('define', help='Define the range by index', parents=[parent_parser])
	def_parser.add_argument('--start', default=None, help='First index in range', type=int, required=True)
	def_parser.add_argument('--stop', default=None, help='Last index in range. Must be larger than start', type=int, required=True)

	args = parent_parser.parse_args()

	# CHECK FOR VALID INPUT
	check = check_input(args)
	if check != 0:
		print('Check input returned non-zero result. Abort.')
		exit()

	outdir = make_output(args)

	if program == 'define':
		start = args.start
		stop = args.stop

	extract_range(args, start, stop, outdir)




if __name__ == "__main__":
	main(sys.argv[1])





