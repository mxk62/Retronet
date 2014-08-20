"""Retrosynthetic worker
"""
from retronet import Chemical


def worker(num, tasks, results, transforms):
    """Performs tasks received from the ventilator.

    Each task is represented by a dictionary

        {'smiles': <smiles>, 'trans_ids': <transform ids> }

    where <smiles> is a string describing a chemical in Daylight's SMILES
    notation and <transform ids> is a list of integers specifying subset of
    transforms which need to be applied.

    Results, all possible reactant sets derived by application of the
    transforms subset, are send back to  the sink also as dictionaries,

        {'results': <reactant_sets>}

    where <reactant_sets> is a list of tuples and each tuple consists of
    strings representing possible reactant set in SMILES notation or is
    None, if none of the transforms from the subset worked.

    Parameters
    ----------
    num : integer
        Worker id number.
    tasks : Queue
        A source of tasks to perform.
    results : Queue
        A sink to put results to.
    transforms : list of Transforms
        A list of available retrosynthetic transforms.
    """
    for msg in iter(tasks.get, 'STOP'):
        smi, ids = msg['smiles'], msg['trans_ids']

        chem = Chemical(smi)

        reactant_sets = []
        for i in ids:
            #print 'Worker-%d: applying transform %d to %s' % (
            #    num, i, chem.smiles)
            reactant_sets.extend(chem.make_retrostep(transforms[i]))

        # Send back the results.
        #print 'Worker-%d: sending back %d results.' % (
        #    num, len(reactant_sets))
        results.put({'results': reactant_sets if reactant_sets else None})
