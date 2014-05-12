from nose.tools import raises
from ..chemical import Chemical
from ..transform import Transform


def test_init_success():
    """Returns True if transform was initialized properly."""

    # Perfectly valid SMILES.
    smiles = 'c1ccccc1'
    chem = Chemical(smiles)
    assert chem.smiles == smiles and chem.mol is not None


@raises(ValueError)
def test_init_failure():
    """Raises ValueError exception due to invalid SMILES."""

    # Invalid SMILES, molecule cannot by kekulized.
    smiles = 'c1cccc1'
    chem = Chemical(smiles)


def test_make_retorstep():
    """Returns True if applying transform gives expected product sets."""

    # Must be a valid transform.
    smarts = '[#6:4][N:2][S:1][#6:3]>>[#7:2][#6:4].[S:1][#6:3]'
    t = Transform(smarts)

    # A valid chemical compatible with transform's retron.
    smiles = 'CSNC(CC1CNSCC1F)C1CSNC(C)C1'
    chem = Chemical(smiles)
    actual_sets = chem.make_retrostep(t)
    expected_sets = {
        frozenset(['CC1CC(C(N)CC2CNSCC2F)CSN1', 'CS']),
        frozenset(['CSNC(CC(CN)CF)C1CSNC(C)C1', 'CSNC(CCC(F)CS)C1CSNC(C)C1']),
        frozenset(['CSNC(CCC(C)N)CC1CNSCC1F', 'CSNC(CC1CNSCC1F)C(C)CS'])
    }
    assert actual_sets == expected_sets

    # A valid chemical NOT compatible with transform's retron.
    smiles = 'c1ccccc1'
    chem = Chemical(smiles)
    actual_sets = chem.make_retrostep(t)
    expected_sets = set()
    actual_sets = chem.make_retrostep(t)
    assert actual_sets == expected_sets
