from rdkit import Chem


class Chemical(object):
    """Represents a chemical compound."""

    def __init__(self, smiles):
        self.smiles = smiles.strip()
        self.mol = Chem.MolFromSmiles(self.smiles)
        if self.mol is None:
            raise ValueError('invalid chemical SMILES: \'%s\'' % self.smiles)

    def __repr__(self):
        return "%s(%r)" % (self.__class__, self.smiles)

    def __str__(self):
        return "%s" % self.smiles

    def has_fragment(self, frag):
        """Returns True if a given fragment is present on chemical.

        Parameters
        ----------
        frag : RDKit query Mol
            A pattern representing the fragment of interest.
        """
        return self.mol.HasSubstructMatch(frag)

    def make_retrostep(self, transform):
        """Returns unique sets of substrate SMILES.

        Parameters
        ----------
        transform : Transform
            A transform which is to applied to chemical.

        Returns
        -------
        smiles : tuples of strings
            A unique set of SMILES representing possible substrates resulting
            from applying the input transform.

            .. warning:: Substrates which could not be sanitized are discarded.
        """
        smis = set()

        products = [self.mol]
        if transform.byproducts:
            products.extend(transform.byproducts)

        reactant_sets = transform.formula.RunReactants(products)
        if reactant_sets:
            for reactants in reactant_sets:
                try:
                    for m in reactants:
                        Chem.SanitizeMol(m)
                except ValueError:
                    continue
                smis.add(frozenset(Chem.MolToSmiles(m) for m in reactants))
        return tuple(tuple(s) for s in smis)