
import re
import collections

file = "../../data/2015-06-08/mykeg.processed"

f = open(file, 'r')

last_A = ''
last_B = ''
last_C = ''

d = collections.defaultdict(dict)
for line in f:
    #print line,
    m = re.match("^A(.*)?", line)
    if m:
        #print m.groups()
        last_A = m.groups()[0].strip()
        last_B = ''
        last_C = ''

    m = re.match("^B(.*)?", line)
    if m:
        #print m.groups()
        last_B = m.groups()[0].strip()
        last_C = ''

    m = re.match("^C(.*)?", line)
    if m:
        #print m.groups()
        last_C = m.groups()[0].strip()

    if last_A and last_B and last_C:
        if last_A in d and last_B in d[last_A]:
            d[last_A][last_B].append(last_C)
        else:
            d[last_A][last_B] = []
            d[last_A][last_B].append(last_C)
f.close()



output_file = "./kegg4R.tsv"
f = open(output_file, 'w')
f.write('A\tB\tC\tpathway\n')

for A in d:
    for B in d[A]:
        for C in d[A][B]:
            m = re.match("(\d+)\s+(.*?)\[PATH:(sce\d+)\]",C)
            if m:
                gr = m.groups()
                f.write('{0}\t{1}\t{2}\tpath:{3}\n'.format(A, B, gr[1], gr[2]))
f.close()
