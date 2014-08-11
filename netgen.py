import argparse
import logging
import multiprocessing
import networkx
import pickle
import time
import zmq
import retronet as rn


def main():
    # Parse command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--depth', type=int, default=1,
                        help='Maximal network depth.')
    parser.add_argument('-p', '--popularity', type=int, default=5,
                        help='Core popularity threshold.')
    parser.add_argument('-t', '--transforms', type=str, default='cores.dat',
                        help='File with available transforms.')
    parser.add_argument('file', type=str,
                        help='File with chemicals SMILES.')
    args = parser.parse_args()

    # Initialize a logger.
    logging.basicConfig(filename='netgen.log',
                        level=logging.INFO,
                        format='%(name)s: %(levelname)s: %(message)s')

    # Acquire all valid, single-retron transforms.
    start = time.clock()
    transforms = get_transforms(args.transforms, popularity=args.popularity)
    eta = time.clock() - start
    logging.info('%d transforms acquired in %f s.' % (len(transforms), eta))

    # Read the SMILES from an external file.
    logging.info('Starting building networks for targets...')
    with open(args.file) as f:
        smiles = [line.strip() for line in f.readlines()]

    # Build the network, either using existing reaction database or by applying
    # retrosynthetic transforms.
    for idx, smi in enumerate(smiles):
        # Time each build individually.
        start = time.clock()
        network = rn.populate(smi, transforms, depth=args.depth)
        eta = time.clock() - start
        logging.info('Finished building for %s in %f s.' % (smi, eta))

        # Gather statistics and write them down to a file.
        stats = rn.count_nodes(network, smi, args.depth)
        line = '\t'.join([str(stats[key]) for key in sorted(stats)])
        with open('counts.dat', 'a') as f:
            f.write('{0}\t{1}\t{2}\n'.format(idx, eta, line))

        # Write graph for debugging purposes
        #networkx.write_gexf(network, 'testgraph.gexf')

        # Discard current network entirely.
        network.clear()
    logging.info('All builds finished.')


def get_transforms(filename, popularity=5):
    """Acquires transforms from an external data file.

    Parameters
    ----------
    filename : string
        Name of the file storing available transforms. File is expected
        to contain a pickled dictionary with keys being transforms encoded
        in Daylight SMARTS format and values representing id numbers of
        reactions transforms applied to.
    popularity : integer
        Popularity threshold, i.e. the minimal number of reaction a
        transform can be applied to. Transforms below the threshold are
        ignored. Defaults to 5.

    Returns
    -------
    transforms : list of Transforms
        A list of transforms.
    """
    with open(filename) as f:
        data = pickle.load(f)
    data = {smarts: rxns
            for smarts, rxns in data.items() if len(rxns) >= popularity}

    transforms = []
    for smarts, rxns in data.items():
        synthons, retrons = smarts.split('>>')
        try:
            t = rn.Transform('{0}>>{1}'.format(retrons, synthons))
        except (RuntimeError, ValueError):
            continue
        if len(t.retrons) == 1:
            transforms.append(t)
    return transforms


if __name__ == '__main__':
    # Create pool of workers to distribute tasks.
    pool_size = multiprocessing.cpu_count()
    pool_size = 1
    for num in range(pool_size):
        w = multiprocessing.Process(target=rn.worker, args=(num,))
        w.start()
    print '%d workers spawned' % pool_size

    context = zmq.Context()

    # Set up a controller.
    controller = context.socket(zmq.PUB)
    controller.bind('tcp://*:5557')

    # Do the stuff.
    main()

    # Terminate workers.
    controller.send_string('DONE')


