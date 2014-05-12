import argparse
import collections
import networkx as nx
import sys
import time
from pymongo import MongoClient
from pymongo.errors import ConnectionFailure
import retronet as rn




# Parse command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument('--ip', type=str, default='127.0.0.1',
                    help='Database server ip.')
parser.add_argument('-d', '--depth', type=int, default=3,
                    help='Maximal network depth.')
parser.add_argument('file', type=str,
                    help='File with chemicals SMILES.')
args = parser.parse_args()

# Connect to the mongo database.
try:
    client = MongoClient(args.ip)
except ConnectionFailure:
    sys.exit(1)
db = client.data

# Acquire all valid, single-retron transforms.
transforms = []
for idx, rec in enumerate(db.retro.find(), start=1):
    try:
        t = rn.Transform(rec['reaction_smarts'].encode('ascii'))
    except (KeyError, ValueError):
        continue
    if len(t.retrons) == 1:
        transforms.append(t)

# Read the SMILES from an external file.
smiles = []
with open(args.file) as f:
    for line in f.readlines():
        smiles.append(line.strip())

# Build the network, either using existing reaction database or by applying
# retrosynthetic transforms.
network = nx.DiGraph()
for idx, smi in enumerate(smiles):
    # Time each build individually.
    start = time.clock()
    network = rn.populate(smi, transforms, depth=args.depth)
    #network = rn.recreate(smi, db, depth=args.depth)
    eta = time.clock() - start

    # Gather statistics and write them down to a file.
    stats = get_stats(network, args.depth, smi)
    line = '\t'.join([str(stats[key]) for key in sorted(stats)])
    with open('occupation.dat', 'a') as f:
        f.write('{0}\t{1}\t{2}\n'.format(idx, eta, line))

    # Write graph for debugging purposes
    #nx.write_gexf(network, 'testgraph.gexf')

    # Discard current network entirely.
    network.clear()
