from nose.tools import raises
from ..transform import Transform


def test_init_success():
    """Returns True if transform was initialized properly."""

    db_id = 4956
    smarts = '[CX4,c:8][C:7]([CX4,c:9])([NH2:2])[C:6]([OH:3])=[O:5].[N:4]>>' \
             '[*:8][C:7]([*:9])=[O:5].[C:6]#[N:4].[O:3].[N:2]'

    # Initializes transform without transform id number.
    t = Transform(smarts)
    assert t.smarts == smarts and t.id is None and t.formula is not None

    # Initializes transform with transform id number.
    t = Transform(smarts, db_id=db_id)
    assert t.smarts == smarts and t.id == db_id and t.formula is not None


@raises(ValueError)
def test_init_failure():
    """Raises an ValueError exception due to improper transform SMARTS."""

    # An invalid transform (missing mapped atom).
    smarts = '[C:7]([CX4,c:9])([NH2:2])[C:6]([OH:3])=[O:5].[N:4]>>' \
             '[*:8][C:7]([*:9])=[O:5].[C:6]#[N:4].[O:3].[N:2]'
    t = Transform(smarts)

    # An invalid transform (empty SMART string).
    smarts = ''
    t = Transform(smarts)
