# xRange
xRange extracts a segment of characters from a DNA sequence file.

## Usage
The user can specify the bases to extract using five different methods. Note that all ranges are inclusive.

Lets say the sequence is `ACTCGTAATGAATAGGTCAGGGGTTCGATTCC`
### Defined coordinates
```
# Extract bases 3-10: TCGTAATG
python xRange.py define -i input.fasta --start 3 --stop 10
```

### Add to position
```
# Extract base 2 and the next 5 bases: CTCGTA
python xRange.py add -i input.fasta --start 2 --add 5
```

### Subtract from position
```
# Extract base 10 and the previous 5 bases: GTAATG
python xRange.py subtract -i input.fasta --stop 10 --sub 5
```

### First x bases
```
# Extract first 10 bases: ACTCGTAATG
python xRange.py first -i input.fasta --num 10
```

### Last x bases
```
#Extract last 10 bases: GTTCGATTCC
python xRange.py last -i input.fasta --num 10
```


### General Usage
By default, xRange extracts a region from the first sequence it sees in the FASTA file. If a file contains multiple sequences, specify one sequence ID using the `--id` flag. To extract the reverse complement of a sequence, use the `-r` flag. Let's say a FASTA file has two sequences: seq1=ACCG and seq2= CCAT
```
# Extract the reverse complement of the first two bases of seq1: GT
python xRange.py first -i input.fasta --num 2 -r

# Extract the first 2 bases of seq2: CC
python xRange.py first -i input.fasta --num 2 --id seq2

```

By default, xRange will create an output direcotry called 'xSequences' located within the current working directory. To change the path to the output directory, use the `-p` or `--output_path` flag. to change the name of this output directory, use the `-o` or `--output_directory` flags.

By default, the name of the exported sequence extract will be `xSeq.fasta`. The prefix of this file can be set using the `-s` or `--savename` flag. The exported sequence will have the same sequence ID as the sequence it was extracted from unless the user wishes to change this using the `--outID` flag. 
