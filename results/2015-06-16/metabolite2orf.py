import csv
import re

file = "./iMM904 model/iMM904_gpr.txt"
output_file = "./iMM904_long.tsv"

f = open(file, 'r')

f_out = open(output_file, 'w')
f_out.write('metabolite\treaction\tec_number\tgene\tside\tdirectionality\n')

reader = csv.reader(f)

for line in reader:
    #print line,

    #parsing reactions
    m = re.search("(-->|<==>)", line[2])
    genes = line[6]
    ec_number = line[4]

    genes = re.sub("[orand\s\(\)]+", " ", genes).strip()
    genes = list(set(genes.split()))

    if not genes:
        continue
    reaction = line[0]

    if m:
        delimiter = m.groups()[0]
        parts = line[2].split(delimiter)
        #print parts
        substrates = parts[0].split("+")
        products = parts[1].split("+")

        #removing trailing spaces
        substrates = [x.strip() for x in substrates]
        products = [x.strip() for x in products]

        #removing stoichiometry
        substrates = [re.sub('\(\d+[\.]?\d*\)\s+', "", x) for x in substrates]
        products = [re.sub('\(\d+[\.]?\d*\)\s+', "", x) for x in products]

        #removing compartments
        substrates = [re.sub('^\[\w+\]\s?:\s?', "", x) for x in substrates]
        products = [re.sub('^\[\w+\]\s?:\s?', "", x) for x in products]

        #removing compartments
        substrates = [re.sub('\[\w+\]', "", x) for x in substrates]
        products = [re.sub('\[\w+\]', "", x) for x in products]

        substrates = list(set(substrates))
        products = list(set(products))

        if sorted(substrates) == sorted(products):
            if not genes:
                continue

        for orf in genes:
            for s in substrates:
                f_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(s, reaction, ec_number, orf, "substrate", delimiter))
            for p in products:
                f_out.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(p, reaction, ec_number, orf, "product", delimiter))


f_out.close()




f.close()
