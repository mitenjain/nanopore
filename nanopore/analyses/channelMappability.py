from nanopore.analyses.abstractAnalysis import AbstractAnalysis
from nanopore.analyses.utils import getFastqDictionary, samIterator
import os
import pysam
import xml.etree.cElementTree as ET
from jobTree.src.bioio import reverseComplement, prettyXml, system
from collections import Counter
import re


class ChannelMappability(AbstractAnalysis):
    """Calculates nanopore channel specific mappability
    """
    def run(self):
        AbstractAnalysis.run(self)
        readSequences = getFastqDictionary(self.readFastqFile)
        nr = re.compile(r"channel_[0-9]+_read_[0-9]+")
        per_channel_read_counts = Counter([int(x.split("_")[1]) for x in readSequences.iterkeys() if re.match(nr, x)])
        sam = pysam.Samfile(self.samFile, "r")
        mapped_read_counts = Counter([int(aR.qname.split("_")[1]) for aR in samIterator(sam) if re.match(nr, aR.qname)])
        assert len(per_channel_read_counts) >= len(mapped_read_counts)
        if len(mapped_read_counts) > 0 and len(per_channel_read_counts) > 0:
            outf = open(os.path.join(self.outputDir, "channel_mappability.tsv"), "w")
            outf.write("Channel\tReadCount\tMappableReadCount\n")
            max_channel = max(513, max(per_channel_read_counts.keys())) #in case there are more than 512 in the future
            for channel in xrange(1, max_channel):
                outf.write("\t".join(map(str, [channel, per_channel_read_counts[channel], mapped_read_counts[channel]])))
                outf.write("\n")
            outf.close()
            system("Rscript nanopore/analyses/channel_plots.R {} {} {}".format(os.path.join(self.outputDir, "channel_mappability.tsv"), os.path.join(self.outputDir, "channel_mappability.pdf"), os.path.join(self.outputDir, "channel_mappability_sorted.png")))
        self.finish()
