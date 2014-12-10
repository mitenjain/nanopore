import xml.etree.cElementTree as ET
import sys

def find_values(f1, f2, f3):
    f1 = ET.parse(f1).getroot()
    f2 = ET.parse(f2).getroot()
    f3 = ET.parse(f3).getroot()

    m1 = float(f1.attrib['avgmismatchesPerReadBase'])
    m2 = float(f2.attrib['avgmismatchesPerReadBase'])
    m3 = float(f3.attrib['avgmismatchesPerReadBase'])

    mismatches = str((m1+m2+m3)/3)

    insert1 = float(f1.attrib['avginsertionsPerReadBase'])
    insert2 = float(f2.attrib['avginsertionsPerReadBase'])
    insert3 = float(f3.attrib['avginsertionsPerReadBase'])

    inserts = str((insert1+insert2+insert3)/3)

    del1 = float(f1.attrib['avgdeletionsPerReadBase'])
    del2 = float(f2.attrib['avgdeletionsPerReadBase'])
    del3 = float(f3.attrib['avgdeletionsPerReadBase'])

    deletes = str((del1+del2+del3)/3)

    i1 = float(f1.attrib['avgidentity'])
    i2 = float(f2.attrib['avgidentity'])
    i3 = float(f3.attrib['avgidentity'])

    identity = str((i1+i2+i3)/3)

    return mismatches, identity, inserts, deletes

results = {}

for l in open(sys.argv[1]):
    l = l.rstrip()
    #analysis = l.split("analysis_")[1].split("/")[0]
    mapper = l.split(".fa_")[1].split("/")[0]
    if mapper not in results:
        results[mapper] = []
    results[mapper].append(l)

for mapper in results:
    f1, f2, f3 = results[mapper]
    results[mapper] = find_values(f1, f2, f3)



with open(sys.argv[2],"w") as f:
    f.write("mapper\tavgMismatch\tavgIdentity\tAvgInsert\tAvgDelete\n")
    for mapper, [mismatches, identity, inserts, delete] in sorted(results.iteritems(), key = lambda x:x[0]):
        if "Realign" in mapper and "Em" not in mapper:
            continue
        else:
            f.write("\t".join([mapper, mismatches, identity, inserts, delete])+"\n")