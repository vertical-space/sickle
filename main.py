import sys

infile = 'example.fastq'
outfile = 'example_trimmed.fastq'

a = open(infile, 'r').read().strip().split('@')
z = open(outfile, 'w')

def trim(Q):
    return len(Q)


for line in a[1:]:
    #print repr(line)
    b = line.strip()
    c = b.split('\n')
    if len(c) == 4:
        Q = c[3]
        d = trim(Q)
        out = '\n'.join(['@'+c[0],c[1][:d+1],c[2],c[3][:d+1]])+'\n'
        z.write(out)


#sys.exit()
