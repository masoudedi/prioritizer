proteininfo = '/home/edi/ngs/04.III/scripts/protein/data/9606.protein.info.v11.5.txt'
proteinalians = '/home/edi/ngs/04.III/scripts/protein/data/9606.protein.aliases.v11.5.txt'
proteinlinkfull = '/home/edi/ngs/04.III/scripts/protein/data/9606.protein.links.full.v11.5.txt'
import json

def getinfo(proteininfo):
    annodict = {}
    with open(proteininfo) as T:
        for line in T:
            line = line.strip('\n').split('\t')
            if len(line) > 1:
                annodict[line[0]] = line[1:]
    return annodict


def getalians(proteinalians):
    annodict = {}
    with open(proteinalians) as T:
        for line in T:
            line = line.strip('\n').split('\t')
            if len(line) > 1:
                annodict[line[1]] = line[0]
    return annodict

def getproteinlink(fulllink, score):
    result = {}
    group = []
    n = 1
    with open(fulllink) as T:
        for line in T:
            line = line.strip('\n').split(' ')
            if len(line) > 1:
                if line[0].startswith('protein1'):
                    result[line[0]] = line[1:]
                else:
                    if int(line[-1]) >= score:
                        if n == 1:
                            Scode = line[0]
                            group.append(line[1:])
                            n += 1
                        else:
                            if line[0] == Scode:
                                group.append(line[1:])
                            else:
                                result[Scode] = group
                                group = []
                                Scode = line[0]
                                group.append(line[1:])
                        # result[line[0]] = line[1:]
        result[Scode] = group
    return result

def geneannotate(fullgene, geneinfo):
    result = {}
    for gene, value in fullgene.items():
        if gene in geneinfo:
            geneanno = geneinfo[gene]
            genename = geneanno[0]
            resultvalue = []
            for valuegene in value:
                if valuegene[0] in geneinfo:
                    valuegeneanno = geneinfo[valuegene[0]]
                    valuegenere = [valuegeneanno[0]]
                    # final = valuegenere #+ valuegene[1:]
                    resultvalue.append(valuegenere[0])
            result[genename] = {"association":resultvalue, "definition":geneanno[2:]}

    return result

geneinfo = getinfo(proteininfo)
fullgene = getproteinlink(proteinlinkfull, 980)
result = geneannotate(fullgene, geneinfo)

# seprator = '|'
# cardinal = '\t'
# output = open('generelationships.txt', 'w')
# for re, val in result.items():
#     newval = [i[0] for i in val[1]]
#     newval = seprator.join(newval)
#     final = [re] + val[0] + [newval]
#     finalre = cardinal.join(final)
#     output.write(finalre)
# output.close()
with open("GeneAssociations.json", "w") as T:
    json.dump(result, T, indent=1)