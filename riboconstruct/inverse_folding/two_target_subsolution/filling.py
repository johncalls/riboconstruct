from ..fillings import filling as fl


class TwoTargetFillingMixin(object):
    def __init__(self):
        self.index = None
        self.unopt = True


class TwoTargetHairpinFilling(fl.HairpinFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, embedded_seq, size, energy):
        fl.HairpinFilling.__init__(
            self, bp_pos_id, bp_id, embedded_seq, size, energy)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetBulgeFilling(fl.BulgeFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, size, energy, pred):
        fl.BulgeFilling.__init__(self, bp_pos_id, bp_id, size, energy, pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetStackFilling(fl.StackFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, energy, pred):
        fl.StackFilling.__init__(self, bp_pos_id, bp_id, energy, pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetInteriorFilling(fl.InteriorFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, embedded_seq, size, energy, pred):
        fl.InteriorFilling.__init__(
            self, bp_pos_id, bp_id, embedded_seq, size, energy, pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetMultiloopFilling(fl.MultiloopFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, free_b_ids_left, num_free_b_left,
                 energy, pred):
        fl.MultiloopFilling.__init__(
            self, bp_pos_id, bp_id, free_b_ids_left, num_free_b_left, energy,
            pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetMultiloopStemFilling(fl.MultiloopStemFilling,
                                    TwoTargetFillingMixin):
    def __init__(self, fixed_bp_id, bp_pos_id, free_b_ids_right,
                 num_free_b_right, energy, stem_filling, pred=None):
        fl.MultiloopStemFilling.__init__(
            self, fixed_bp_id, bp_pos_id, free_b_ids_right, num_free_b_right,
            energy, stem_filling, pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetExteriorFilling(fl.ExteriorFilling, TwoTargetFillingMixin):
    def __init__(self, fixed_bp_id, bp_pos_id, free_b_ids_right,
                 num_free_b_right, energy, stem_filling, pred=None):
        fl.ExteriorFilling.__init__(
            self, fixed_bp_id, bp_pos_id, free_b_ids_right, num_free_b_right,
            energy, stem_filling, pred)
        TwoTargetFillingMixin.__init__(self)


class TwoTargetFinalFilling(fl.FinalFilling, TwoTargetFillingMixin):
    def __init__(self, bp_pos_id, bp_id, free_b_id_left, num_free_b_left,
                 energy, pred):
        fl.FinalFilling.__init__(
            self, bp_pos_id, bp_id, free_b_id_left, num_free_b_left, energy,
            pred)
        TwoTargetFillingMixin.__init__(self)
