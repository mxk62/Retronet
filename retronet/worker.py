"""

"""
from retronet import Chemical


def worker(num, tasks, results, transforms):
    """Performs retrosynthetic steps.

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
