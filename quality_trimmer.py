"""
quality_trimmer.py

usage:
python quality_trimmer.py infile x y z

where:
infile = input file to be trimmed. Must be fastq format (preferably .fastq suffix)
x = minimum average quality across window, (Phred score), e.g. 57
y = minimum length threshold after trimming, e.g. 30
z = window size, e.g. 10

example:
python quality_trimmer.py example.fastq 57 30 10

"""

from sys import argv

infile = argv[1]	# input file. Must be fastq format, preferably .fastq
outfile = argv[1].rstrip('.fastq')+'_trimmed.fastq'

qual = int(argv[2])	# minimum average quality across window, (Phred score), e.g. 57
qual_min = qual + 33

len_min = int(argv[3])	# minimum length threshold after trimming, e.g. 30
win = int(argv[4])	# window size, e.g. 10

a = open(infile, 'r').read().strip().split('@')
z = open(outfile, 'w')

discard = 0
trimmed = 0
intact = 0

def trim(Q):
    global qual_min, len_min, win
    endpoint = len(Q)-win
    out = ''
    qualscores = [ord(chr) for chr in Q]
    x = 0
    while True:
        average = sum(qualscores[x:x+win])/float(win)
        if average < qual_min or x > endpoint:
            break
        x += 1
    if x+win > len_min:
        return x+win
    else:
        return None


for line in a[1:]:
    b = line.strip()
    c = b.split('\n')
    if len(c) == 4:
        Q = c[3]
        d = trim(Q)
        if d != None:
            out = '\n'.join(['@'+c[0],c[1][:d+1],c[2],c[3][:d+1]])+'\n'
            if d < len(Q):
                trimmed += 1
            else:
                intact += 1
            z.write(out)
        else:
            discard += 1

print 'intact reads:', intact
print 'trimmed reads:', trimmed
print 'discarded reads:', discard

assert discard+trimmed+intact == len(a)-1

