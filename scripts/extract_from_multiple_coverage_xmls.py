import xml.etree.cElementTree as ET
import sys
f1 = ET.parse(sys.argv[1]).getroot()
f2 = ET.parse(sys.argv[2]).getroot()
f3 = ET.parse(sys.argv[3]).getroot()

lengths = list()
identity = list()
coverage = list()
insertions = list()
deletions = list()
mismatches = list()


outf = open(sys.argv[4],"w")

for c in f1.getchildren():
    lengths.append(c.attrib["readLength"])
    identity.append(c.attrib["identity"])
    insertions.append(c.attrib["insertionsPerReadBase"])
    deletions.append(c.attrib["deletionsPerReadBase"])
    mismatches.append(c.attrib["mismatchesPerReadBase"])

for c in f2.getchildren():
    lengths.append(c.attrib["readLength"])
    identity.append(c.attrib["identity"])
    insertions.append(c.attrib["insertionsPerReadBase"])
    deletions.append(c.attrib["deletionsPerReadBase"])
    mismatches.append(c.attrib["mismatchesPerReadBase"])

for c in f3.getchildren():
    lengths.append(c.attrib["readLength"])
    identity.append(c.attrib["identity"])
    insertions.append(c.attrib["insertionsPerReadBase"])
    deletions.append(c.attrib["deletionsPerReadBase"])
    mismatches.append(c.attrib["mismatchesPerReadBase"])

outf.write("length ")
outf.write(" ".join(lengths)+"\n")
outf.write("identity ")
outf.write(" ".join(identity)+"\n")
outf.write("insertions ")
outf.write(" ".join(insertions)+"\n")
outf.write("deletions ")
outf.write(" ".join(deletions)+"\n")
outf.write("mismatches ")
outf.write(" ".join(mismatches)+"\n")
outf.close()