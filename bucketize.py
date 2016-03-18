import sys

source = sys.argv[1]

d = dict()

f = open("buckets.list")
for r in f:
	r = r.rstrip()
	seq,b = r.split("\t")
	if seq in d:
		d[seq].append(b)
	else:
		d[seq] = [b]

f.close()

f = open(source)
buckets = []
for r in f:
	seq = r.rstrip()
	if seq[0]==">":
		for b in buckets:
			b.close()
		buckets = list(map(lambda x: open(x+".fasta","a"), d.get(seq[1:],[])))
	for b in buckets:
		b.write(r)
for b in buckets:
	b.close()
f.close()
		
