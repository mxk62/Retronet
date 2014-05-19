import chemical as chem


class Reaction(object):
    """Represents a reaction."""

    def __init__(self, smiles):
        """Initialize entry."""
        self.smiles = smiles.strip()
        self.react_smis, self.prod_smis = [s.split('.')
                                           for s in self.smiles.split('>>')]
        try:
            self.reactants = [chem.Chemical(smi) for smi in self.react_smis]
            self.products = [chem.Chemical(smi) for smi in self.prod_smis]
        except ValueError:
            raise ValueError('invalid substrates and/or products')

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles
