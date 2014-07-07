from rdkit.Chem import AllChem


class Transform(object):
    """Represents retrosynthetic transform."""

    def __init__(self, smarts, db_id=None):
        self.id = db_id
        self.smarts = smarts
        self.retrons, self.synthons = [s.split('.')
                                       for s in smarts.split('>>')]
        self.formula = AllChem.ReactionFromSmarts(smarts)
        if self.formula is None:
            raise ValueError('invalid transform SMARTS')

    def __repr__(self):
        return "%s(%r, db_id=%r)" % (self.__class__, self.smarts, self.id)

    def __str__(self):
        return "%s" % self.smarts
