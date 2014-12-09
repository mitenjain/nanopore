import sys, os

readTypes = ["2D","complement","template"]
files = []
for f in os.listdir(sys.argv[1]):
    if f.split("_")[0] in readTypes:
        readType = f.split("_")[0]
    if "blast_report.txt" in f and "Fail" not in f:
        files.append([os.path.join(sys.argv[1], f), readType])
files = sorted(files, key=lambda x: x[1])

with open(sys.argv[2],"w") as f:
    for infile, readType in files:
        readType = readType.title()
        start = """
\clearpage \n
\\begin{{longtable}}{{|p{{14cm}}|p{{1.1cm}}|}} \hline \n
\multicolumn{{2}}{{|c|}}{{{} Unmapped Reads BLAST Hits}}  \\\\ \n
Sequence Name & Counts \\\\ \hline \n""".format(readType)

        end = """
\hline \n
\caption{{Table of BLAST hits for {0} reads unmapped by any mapper.}} \n
\label{{{0}BlastTable}} \n
\end{{longtable}} \n""".format(readType)

        data = [x.split("\t")[-2:] for x in open(infile)][1:]
        tmp = []
        for x in data:
            x[0].replace("_","\_")
            tmp.append(x)

        f.write(start)
        for x in tmp:
            f.write(" & ".join(x)+" \\\\ \n")

        f.write(end)
