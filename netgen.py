import argparse
import logging
import networkx
import pymongo
import time
import retronet as rn


def main():
    # Parse command line arguments.
    parser = argparse.ArgumentParser()
    parser.add_argument('-C', '--collection', type=str, default='retros',
                        help='Database collection holding transforms, '
                             'defaults to \'retros\'.')
    parser.add_argument('-D', '--database', type=str, default='nerc',
                        help='Database to connect to, defaults to \'nerc\'.')
    parser.add_argument('-p', '--popularity', type=int, default=5,
                        help='Transform popularity threshold ')
    parser.add_argument('-d', '--depth', type=int, default=1,
                        help='Maximal network depth.')
    parser.add_argument('file', type=str,
                        help='File with chemicals SMILES.')
    args = parser.parse_args()

    # Initialize a logger.
    logging.basicConfig(filename='netgen.log',
                        level=logging.INFO,
                        format='%(name)s: %(levelname)s: %(message)s')

    # Acquire all valid transforms.
    start = time.clock()
    client = pymongo.MongoClient()
    db = client[args.database]
    transforms = get_transforms(db, args.collection,
                                popularity=args.popularity)
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
        start = time.time()
        network = rn.populate(smi, transforms, depth=args.depth)
        eta = time.time() - start
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


def get_transforms(db, collection, popularity=5):
    """Acquires transforms from an external data file.

    Parameters
    ----------
    db : mongodb
        Name of a mongodb.
    collection : mongodb collection
        Collection containing retrosynthetic transforms.
    popularity : integer
        Popularity threshold, i.e. the minimal number of reaction a
        transform can be applied to. Transforms below the threshold are
        ignored. Defaults to 5.

    Returns
    -------
    transforms : list of Transforms
        A list of transforms.
    """
    transforms = []
    for rec in db[collection].find({'expert': True}):
        byproducts = None
        if 'product_smiles' in rec:
            byproducts = [s.encode('ascii') for s in rec['product_smiles']]

        if 'popularity' in rec and rec['populairty'] < popularity:
            continue

        try:
            t = rn.Transform(rec['reaction_smarts'].encode('ascii'),
                             db_id=rec['_id'], byproducts=byproducts)
        except (KeyError, ValueError):
            continue

        # Some transforms has to many byproducts listed (e.g. 87), filter
        # them out.
        nb = len(t.byproducts) if t.byproducts is not None else 0
        if nb == len(t.retrons) - 1:
            transforms.append(t)
    return transforms


if __name__ == '__main__':
    main()
