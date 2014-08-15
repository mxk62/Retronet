from retronet import Chemical


def worker(num, in_queue, out_queue, transforms):
    """Performs retrosynthetic steps.

    Parameters
    ----------
    num : integer
        Worker id number.
    in_queue : Queue
        Source of steps to perform
    out_queue : Queue
        Sink to put results to.
    transforms : list of Transforms
        List of available transforms.
    """
    for msg in iter(in_queue.get, 'STOP'):
        smiles, ids = msg['smiles'], msg['trans_ids']

        chem = Chemical(smiles)

        reactant_sets = []
        for i in ids:
            #print 'Worker-%d: applying transform %d to %s' % (
            #    num, i, chem.smiles)
            reactant_sets.extend(chem.make_retrostep(transforms[i]))

        # Send back the results.
        #print 'Worker-%d: sending back %d results.' % (
        #    num, len(reactant_sets))
        out_queue.put({'results': reactant_sets if reactant_sets else 'NONE'})
