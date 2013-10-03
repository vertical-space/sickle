quality_trimmer.py is a command line application for quality trimming of fastq files. 

usage:
python quality_trimmer.py infile x y z

where:
infile = input file to be trimmed. Must be fastq format (preferably with .fastq suffix)
x = minimum average quality across sliding window, (Phred score), e.g. 57
y = minimum length after trimming to keep read, e.g. 30
z = window size, e.g. 10

example session (at the command line prompt):
python quality_trimmer.py example.fastq 57 30 10

# note: failed sequencess are excluded from the final output
