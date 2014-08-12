from rdkit import Chem
from rdkit.Chem import AllChem


class Transform(object):
    """Represents retrosynthetic transform."""

    def __init__(self, smarts, db_id=None, byproducts=None, explicit_hs=False):
        self.id = db_id

        # List of acceptable byproducts that can be supplied implicitly.
        self.byproducts = []
        if byproducts is not None:
            self.byproducts.extend([Chem.MolFromSmiles(s) for s in byproducts])

        # Actual chemical reaction (an instance of RDKit's ChemicalReaction)
        self.formula = AllChem.ReactionFromSmarts(smarts)
        if self.formula is None:
            raise ValueError('invalid transform SMARTS')

        self.retrons, self.synthons = [s.split('.')
                                       for s in smarts.split('>>')]
        self.smarts = smarts
        self.explicit_hs = explicit_hs

    def __repr__(self):
        return "%s(%r, db_id=%r)" % (self.__class__, self.smarts, self.id)

    def __str__(self):
        return "%s" % self.smarts
