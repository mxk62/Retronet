from rdkit import Chem


class Chemical(object):
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smiles = smiles.strip()
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError('invalid SMILES')

    def make_retrostep(self, transform):
        """Returns unique sets of substrate SMILES.

        Function returns unique sets of SMILES representing possible
        substrates obtained by applying a given transform to the product.
        """
        smis = set()
        product_sets = transform.formula.RunReactants((self.mol,))
        if product_sets:
            for products in product_sets:
                try:
                    for m in products:
                        Chem.SanitizeMol(m)
                except ValueError:
                    continue
                smis.add(frozenset([Chem.MolToSmiles(m) for m in products]))
        return smis