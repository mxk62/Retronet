from rdkit import Chem


class Chemical(object):
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smiles = smiles.strip()
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError('invalid SMILES')

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles

    def has_fragment(self, frag):
        """Returns True if a given fragment is present on chemical.

        Parameters
        ----------
        frag : Pattern
            A pattern representing the fragment of interest.
        """
        return self.mol.HasSubstructMatch(frag.templates[0])

    def make_retrostep(self, transform):
        """Returns unique sets of substrate SMILES.

        Parameters
        ----------
        transform : Transform
            A transform which is to applied to chemical.

        Returns
        -------
        smiles : set
            Unique set of SMILES representing possible substrates resulting
            from applying the input transform.

            .. warning:: Substrates which could not be sanitized are discarded.
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