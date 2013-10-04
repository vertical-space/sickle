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

def trim(Q, qual_min, len_min, win):
    """
    This is the core function. 

    It takes as input the quality string from each sequence in the fastq file.

    it implements a sliding window (of user-specified length), which reads from left to right 
    and tests whether the average score across the window is above the user-specified cutoff.

    It returns the index corresponding to the base that caused the average score to drop below the threshold 
    (or the last base if this did not occur), which indicates the position where the read needs to 
    be cut, or returns None if the length after trimming is below the user-specified length cutoff.

    """
    endpoint = len(Q)-win 
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


def main():
    infile = argv[1]	# input file. Must be fastq format, preferably .fastq
    outfile = argv[1].rstrip('.fastq')+'_trimmed.fastq' #construct outfile name
    
    qual = int(argv[2])	# minimum average quality across window, (Phred score), e.g. 57
    qual_min = qual + 33	# fastq scores are encoded as Phred+33
    
    len_min = int(argv[3])	# minimum length threshold after trimming, e.g. 30
    win = int(argv[4])	# window size, e.g. 10
    
    a = open(infile, 'r').read().strip().split('@')
    z = open(outfile, 'w')
    
    # initialize counters to track progress
    discard = 0
    trimmed = 0
    intact = 0
    
    # process the file
    for line in a[1:]:	# item 0 is an empty string so ignore it
        b = line.strip()
        c = b.split('\n')
        if len(c) == 4:
            Q = c[3]	# this is the quality string
            d = trim(Q, qual_min, len_min, win)
            if d != None:
                # construct the printline in fastq format:
                out = '\n'.join(['@'+c[0],c[1][:d+1],c[2],c[3][:d+1]])+'\n'	
                if d < len(Q):
                    trimmed += 1
                else:
                    intact += 1
                # write to file:
                z.write(out)
            else:
                # note: failed sequencess are excluded from the final output
                discard += 1
    
    print 'intact reads:', intact
    print 'trimmed reads:', trimmed
    print 'discarded reads:', discard
    
    assert discard+trimmed+intact == len(a)-1

if __name__ == "__main__":
    main()