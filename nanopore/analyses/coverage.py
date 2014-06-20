from nanopore.analyses.abstractAnalysis import AbstractAnalysis
import os

class Coverage(AbstractAnalysis):
    """This is just my first test analysis target
    """
    def run(self):
        outputFile = os.path.join(self.outputDir, "coverage.tsv")
        open(outputFile, 'w').close()