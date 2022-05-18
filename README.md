# xRange
Extract a segment of characters from a DNA sequence file
A script that provides several options for extracting a sequence from a FASTA file


## Usage
There are three ways to specify the range of sequence to extract: defined coordinates, add from position, and subtract from position
### Defined coordinates
```
# Use defined coordinates to extract first 10 bases
python xRange.py define -f file.fasta --start 1 --stop 10

# Output sequence to file fist.fasta
python xRange.py define -f file.fasta --start 1 --stop 10 -s first

# Set read ID
python xRange.py define -f file.fasta --start 1 --stop 10 -s first --outID myID
```

### Add to position


### Subtract from position
```
# Extract positions 1 through 4
python xRange.py subtract -f tmp.fa --stop 4 --sub 2
```
