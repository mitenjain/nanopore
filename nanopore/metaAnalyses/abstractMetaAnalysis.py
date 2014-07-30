from jobTree.scriptTree.target import Target

class AbstractMetaAnalysis(Target):
    """Base class to for meta-analysis targets. Inherit this class to create a meta-analysis.
    """
    def __init__(self, outputDir, experiments):
        Target.__init__(self)
        self.experiments = experiments
        self.outputDir = outputDir
