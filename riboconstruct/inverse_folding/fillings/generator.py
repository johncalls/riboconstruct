from ... import rna


class FillingsGenerator(object):
    def __init__(self, subsolution):
        self.subsolution = subsolution

    def __getattr__(self, attr):
        """Everything else is delegated to the object"""
        return getattr(self.subsolution, attr)

    def evaluate(self):
        for bp_id in xrange(rna.BasepairId.count):
            self.generate_fillings(bp_id)
        self.subsolution.delete_predecessor()
