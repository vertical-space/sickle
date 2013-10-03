qual = 57	# variable (Phred score)

qual_min = qual + 33
len_min = 30
win = int(10)

infile = 'example.fastq'
outfile = 'example_trimmed.fastq'

a = open(infile, 'r').read().strip().split('@')
z = open(outfile, 'w')

discard = 0
trimmed = 0
intact = 0

def trim1(Q):
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

def trim2(Q):
    global qual_min, len_min, win
    endpoint = len(Q)-win
    out = ''
    qualscores = [ord(chr) for chr in Q]
    x = 0
    while True:
        #print x,endpoint
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
        d = trim2(Q)
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
