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
    trans_graph.add_nodes_from([p for p in patterns])
else:
    trans_graph = rn.create_depgraph(patterns)

# Having sequences as node attributes make write_gexf() fail, thus we are
# making a copy of the graph with lists of duplicates stripped.
g = trans_graph.copy()
for v in g:
    try:
        del g.node[v]['duplicates']
    except KeyError:
        continue

# Write down the dependency graph.
networkx.write_gexf(g, 'depgraph.gexf')

# Write down example chemical and its descendants.
chem_env = networkx.bfs_tree(g, 'c1ccccc1')
networkx.write_gexf(chem_env, 'benzene.gexf')

# Associate transforms with graph patterns
for patt in trans_graph.nodes():
    trans_graph.node[patt]['transforms'] = cores[patt.smiles]
    try:
        for dup in trans_graph.node[patt]['duplicates']:
            trans_graph.node[patt]['transforms'].extend(cores[dup])
    except KeyError:
        continue

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
