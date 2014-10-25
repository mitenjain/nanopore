"""Lib for generating latex tables.
"""

def formatFloat(string, decimals=3):
    f = float(string)
    if f == 2147483647:
        return "NaN"
    return (("%." + str(decimals) + "f") % f)[1:]

def writeDocumentPreliminaries(fileHandle):
    fileHandle.write("\\documentclass[10pt]{article}\n")
    
    fileHandle.write("\\pagenumbering{arabic}\n")
    fileHandle.write("\\pagestyle{plain}\n")
    
    fileHandle.write("\\usepackage{epsfig}\n")
    fileHandle.write("\\usepackage{url}\n")
    fileHandle.write("\\usepackage{rotating}\n")
    
    fileHandle.write("\\usepackage{multirow}\n")
    
    fileHandle.write("\\usepackage{color}\n")
    fileHandle.write("\\usepackage[table]{xcolor}\n")
    
    fileHandle.write("\\setlength{\\evensidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\oddsidemargin}{0in}\n")
    fileHandle.write("\\setlength{\\marginparwidth}{1in}\n")
    fileHandle.write("\\setlength{\\textwidth}{6.5in}\n")
    fileHandle.write("\\setlength{\\topmargin}{-0.5in}\n")
    fileHandle.write("\\setlength{\\textheight}{9in}\n")

    fileHandle.write("\\begin{document}\n")
    
def writeDocumentEnd(fileHandle):
    fileHandle.write("\\end{document}\n")
 
def writePreliminaries(columnNumber, fileHandle):
    fileHandle.write("\\begin{sidewaystable}\n\\centering\n")
    fileHandle.write("\\begin{tabular}{" + ("c"*columnNumber) + "}\n")

def writeEnd(fileHandle, tableLabel, caption):
    fileHandle.write("\\end{tabular}\n")
    fileHandle.write("\caption{%s}\n" % caption)
    fileHandle.write("\label{%s}\n" % tableLabel)
    fileHandle.write("\end{sidewaystable}\n\n")

def writeLine(columnNumber, rowNumber, entries, fileHandle, trailingLines=1):
    updatedEntries = []
    for name, x1, x2, y1, y2 in entries:
        columnNumber = y2 - y1 + 1
        updatedEntries.append((y1, x1, x2, name, columnNumber, y2 - y1 == 0))
        while y2 - y1 > 0: #Make multiple new entries:
            y1 += 1
            updatedEntries.append((y1, x1, x2, "", columnNumber, y2 - y1 == 0))
    updatedEntries.sort()
    #fileHandle.write("\hline\n")
    start = True
    currentRow = 0
    row = []
    for y1, x1, x2, name, columnNumber, cLine in updatedEntries:
        if y1 != currentRow:
            fileHandle.write(" \\\\ %s\n" % " ".join([ "\\cline{%i-%i}" % (x3+1, x4+1) for x3, x4 in row ]))
            currentRow = y1
            row = []
        else:
            if not start:
                fileHandle.write(" & ")
            else:
                start = False
        if cLine:
            row.append((x1, x2))
        fileHandle.write("\multicolumn{%i}{c}{\multirow{%i}{*}{%s}}" % (x2-x1+1, columnNumber, name))
    fileHandle.write(" \\\\\n")
    for i in xrange(trailingLines):
        fileHandle.write("\hline\n")

def writeRow(entries, fileHandle):
    fileHandle.write("%s \\\\\n" % " & ".join(entries))
    
def writeFigure(fileHandle, imageFile, caption, label, width=10):
    fileHandle.write("\\clearpage\n")
    fileHandle.write("\\begin{figure}[h!]\n\\begin{center}\n\\includegraphics[width=%scm]{%s}\n\\caption{%s}\n\\label{%s}\n\\end{center}\n\\end{figure}\n\n" % (width, imageFile, caption, label))
