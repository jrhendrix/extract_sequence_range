# xRange
A script that provides several options for extracting a range of characters from a FASTA file


## Usage
### Specify range with defined coordinates
```
# Use defined coordinates to extract first 10 bases
python xRange.py define -f file.fasta --start 1 --stop 10

# Output sequence to file fist.fasta
python xRange.py define -f file.fasta --start 1 --stop 10 -s first

# Set read ID
python xRange.py define -f file.fasta --start 1 --stop 10 -s first --outID myID
```
