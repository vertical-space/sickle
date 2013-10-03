qual = 47	# variable (Phred score)

qual_min = qual + 33
len_min = 20

infile = 'example.fastq'
outfile = 'example_trimmed.fastq'

a = open(infile, 'r').read().strip().split('@')
z = open(outfile, 'w')

discard = 0

def trim(Q):
    global qual_min, len_min
    out = ''
    for chr in Q:
        if ord(chr) < qual_min:
            break
        out += chr
    if len(out) > len_min:
        #print out
        return len(Q)
    else:
        return None


for line in a[1:]:
    #print repr(line)
    b = line.strip()
    c = b.split('\n')
    if len(c) == 4:
        Q = c[3]
        d = trim(Q)
        if d != None:
            out = '\n'.join(['@'+c[0],c[1][:d+1],c[2],c[3][:d+1]])+'\n'
            print out,
            z.write(out)
        else:
            discard += 1

print 'discarded reads:', discard
