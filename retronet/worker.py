import time
from retronet import Chemical


def worker(num, connection, transforms):
    """Performs retrosynthetic steps.

    Parameters
    ----------
    num : integer
        Worker id number.
    connection : Connection
        Connection to ventilator.
    transforms : list of Transforms
        List of available transforms.
    """
    while connection.poll():
        msg = connection.recv()

        print 'Worker-%d: received \'%s\'' % (num, msg)
        time.sleep(1)

        smiles, ids = msg['smiles'], msg['trans_ids']

        chem = Chemical(smiles)

        reactant_sets = []
        for i in ids:
            reactant_sets.extend(chem.make_retrostep(transforms[i]))

        # Send back the results.
        connection.send(
            {'results': reactant_sets if reactant_sets else 'NONE'}
        )

        print 'Worker-%d: sent back %d results.' % (num, len(reactant_sets))
        time.sleep(1)

