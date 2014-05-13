import argparse
import sys
import pickle
import time
import networkx
import retronet as rn


# Parse command line arguments.
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--depth', type=int, default=3,
                    help='Maximal network depth.')
parser.add_argument('-t', '--transforms', type=str, default='cores.dat',
                    help='File with available transforms.')
parser.add_argument('-p', '--popularity', type=int, default=5,
                    help='Core popularity threshold.')
parser.add_argument('file', type=str,
                    help='File with chemicals SMILES.')
args = parser.parse_args()

# Acquire all valid, single-retron transforms.
with open(args.transforms) as f:
    data = pickle.load(f)
transforms = []
for smarts, rxns in data.items():
    if len(rxns) >= args.popularity:
        synthons, retrons = smarts.split('>>')
        try:
            t = rn.Transform('{0}>>{1}'.format(retrons, synthons))
        except ValueError:
            continue
        if len(t.retrons) == 1:
            transforms.append(t)

cores = {}
patterns = set()
for t in transforms:
    try:
        patt = rn.Pattern(t.retrons[0])
    except ValueError:
        continue
    patterns.add(patt)
    cores.setdefault(patt.smiles, []).append(t)
if False:
    trans_graph = networkx.DiGraph()
    for patt in patterns:
        trans_graph.add_node(patt.smiles)
else:
    trans_graph = rn.create_depgraph(patterns)

# Write down the dependency graph.
networkx.write_gexf(trans_graph, 'depgraph.gexf')

# Write down example chemical and its descendants.
chem_env = networkx.bfs_tree(trans_graph, 'c1ccccc1')
networkx.write_gexf(chem_env, 'benzene.gexf')

print 'Graph size: ', len(trans_graph.nodes()), len(trans_graph.edges())

# Associate transforms with graph patterns
for smi in trans_graph.nodes():
    trans_graph.node[smi]['transforms'] = cores[smi]
    try:
        for dup in trans_graph.node[smi]['duplicates']:
            trans_graph.nodes[smi]['transforms'].extend(cores[dup])
    except KeyError:
        continue

print '# transforms:', sum(len(trans_graph.node[v]['transforms'])
                           for v in trans_graph.nodes())

# Read the SMILES from an external file.
with open(args.file) as f:
    smiles = [line.strip() for line in f.readlines()]

# Build the network, either using existing reaction database or by applying
# retrosynthetic transforms.
for idx, smi in enumerate(smiles):
    # Time each build individually.
    start = time.clock()
    network = rn.populate(smi, trans_graph, depth=args.depth)
    #network = rn.recreate(smi, db, depth=args.depth)
    eta = time.clock() - start

    # Gather statistics and write them down to a file.
    stats = rn.count_nodes(network, smi, args.depth)
    line = '\t'.join([str(stats[key]) for key in sorted(stats)])
    with open('occupation.dat', 'a') as f:
        f.write('{0}\t{1}\t{2}\n'.format(idx, eta, line))

    # Write graph for debugging purposes
    networkx.write_gexf(network, 'testgraph.gexf')

    # Discard current network entirely.
    network.clear()
