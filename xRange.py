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

	# INDEX MODIFIER and FIX OUT OF BOUNDS CASES
	## Note: Input will be position. Translate to index
	if start != None:
		start = start - 1
		if start < 0:
			start = 0
	if stop != None:
		stop = stop - 1
		if stop < 0:
			stop = 0

	if stop is None or start is None:
		return start, stop

	# CHECK THAT 'STOP' COMES AFTER 'START'
	if stop < start: # flip indecies
		print('Start index is larger than stop index. Inverting...\n')
		iStart = start
		start = stop
		stop = iStart

	return start, stop




def extract_range(args, start, stop, flag, outdir):
	'''Find a range of indexes and write bases to file'''

	start, stop = modify_index(start, stop)

	outname = ''.join((args.savename, '.fasta'))
	outfile = '/'.join((outdir.path, outname))
	f2 = open(outfile, 'w')

	# LOOP THROUGH SEQUENCE FILE
	records = []
	f1 = open(args.input_fasta, 'r')
	for record in SeqIO.parse(f1, 'fasta'):
		# Restrict extraction to a specific node
		if args.id != None and args.id != record.id:
			continue

		max_index = len(record.seq)

		# INDEX MODIFIER
		if stop is None:
			stop = max_index
		elif stop > max_index:
			stop = max_index

		# ADJUST START FOR LAST X METHOD
		if flag == 'last':
			start = stop - args.num
			if start < 0:
				start = 0


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
	if not args.input_fasta:
		print("ERROR: Could not find file")
		return 1
	try:
		f = File(args.input_fasta)
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
	d = Dir(args.output_path)
	outpath = '/'.join((args.output_path, args.output_directory))
	if os.path.isdir(outpath):
		outdir = Dir(outpath)
	else:
		os.makedirs(outpath)
		outdir = Dir(outpath)

	return outdir

def define(args):
	start = args.start
	stop  = args.stop
	return start, stop, 'define'

def add(args):
	start = args.start
	stop = start + args.add
	return start, stop, 'add'

def sub(args):
	stop = args.stop
	start = stop - args.sub
	return start, stop, 'sub'

def first(args):
	start = 1
	stop = args.num
	return start, stop, 'first'

def last(args):
	start = None
	stop = None
	return start, stop, 'last'

def main(args):

	command = 'Command: %s' % ' '.join(sys.argv)
	print(f'Running: ', command)

	# CHECK FOR VALID INPUT
	check = check_input(args)
	if check != 0:
		print('Check input returned non-zero result. Abort.')
		exit()

	# SET UP OUTPUT
	outdir = make_output(args)

	# GET START AND STOP
	start, stop, flag = args.func(args)

	# EXTRACT SEQUENCE
	extract_range(args, start, stop, flag, outdir)



if __name__ == "__main__":
	cwd = os.getcwd()

	# PARSER : ROOT
	parent_parser = argparse.ArgumentParser(prog='xRange')
	subparsers = parent_parser.add_subparsers(help='available actions')
	subparsers.required = True

	# DEFINE SUBPARSERS
	parser_def = subparsers.add_parser('define', help='Define the range by index')
	parser_def.set_defaults(func=define)

	parser_add = subparsers.add_parser('add', help='Define the range by first index')
	parser_add.set_defaults(func=add)

	parser_sub = subparsers.add_parser('subtract', help='Define the range by last index')
	parser_sub.set_defaults(func=sub)

	parser_first = subparsers.add_parser('first', help='Extract the first X bases')
	parser_first.set_defaults(func=first)

	parser_last = subparsers.add_parser('last', help='Extract the last X bases')
	parser_last.set_defaults(func=last)

	# PARSER : DEFINED
	
	parser_def.add_argument('-i', '--input_fasta', help="FASTA file containing sequences")
	parser_def.add_argument('--id', default=None, help='Match sequence ID')
	parser_def.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parser_def.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parser_def.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_def.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parser_def.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	parser_def.add_argument('--start', default=None, help='First index in range', type=int, required=True)
	parser_def.add_argument('--stop', default=None, help='Last index in range. Must be larger than start', type=int, required=True)

	
	# PARSER : ADD
	parser_add.add_argument('--add', default=None, help='Length of extract', type=int, required=True)
	parser_add.add_argument('-i', '--input_fasta', help="FASTA file containing sequences")
	parser_add.add_argument('--id', default=None, help='Match sequence ID')
	parser_add.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parser_add.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parser_add.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_add.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parser_add.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	parser_add.add_argument('--start', default=None, help='First index in range', type=int, required=True)

	
	# PARSER : SUBTRACT
	parser_sub.add_argument('-i', '--input_fasta', help="FASTA file containing sequences")
	parser_sub.add_argument('--id', default=None, help='Match sequence ID')
	parser_sub.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parser_sub.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parser_sub.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_sub.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parser_sub.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	parser_sub.add_argument('--stop', default=None, help='Last index in range', type=int, required=True)
	parser_sub.add_argument('--sub', default=None, help='Length of extract', type=int, required=True)


	# PARSER : FIRST X
	parser_first.add_argument('-i', '--input_fasta', help="FASTA file containing sequences")
	parser_first.add_argument('--id', default=None, help='Match sequence ID')
	parser_first.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parser_first.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parser_first.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_first.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parser_first.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	parser_first.add_argument('--num', default=None, help='Number of bases to extract', type=int, required=True)
	
	# PARSER : LAST X
	parser_last.add_argument('-i', '--input_fasta', help="FASTA file containing sequences")
	parser_last.add_argument('--id', default=None, help='Match sequence ID')
	parser_last.add_argument('-o', '--output_directory', default='xSequences', help='Prefix of output directory', type=str)
	parser_last.add_argument('--outID', default=None, help='ID for output sequence', type=str)
	parser_last.add_argument('-p', '--output_path', default=cwd, help='Path to output', type=str)
	parser_last.add_argument('-r', '--reverse_complement', default=False, action='store_true', help='Gets the reverse compliment')
	parser_last.add_argument('-s', '--savename', default='xSeq', help='Prefix for output file')
	parser_last.add_argument('--num', default=None, help='Number of bases to extract', type=int, required=True)

	args = parent_parser.parse_args()

	main(args)





