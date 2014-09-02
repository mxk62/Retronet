import argparse
import itertools
import jsonschema
import pickle
from pymongo import MongoClient
from rdkit.Chem import AllChem
from schema import SCHEMA


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--database', type=str, default='test',
                        help='Database to put records into.')
    parser.add_argument('-c', '--collection', type=str, default='auto',
                        help='Collection to put records into.')
    parser.add_argument('filename', type=str, help='File with transforms.')
    args = parser.parse_args()

    rec_count = itertools.count(start=1)
    records = []
    transforms = get_transforms(args.filename)
    for rxn_smarts, rxns in transforms.items():
        retro_smarts = reverse_smarts(rxn_smarts)
        try:
            rxn = AllChem.ReactionFromSmarts(retro_smarts)
        except RuntimeError:
            continue

        rec = make_record(retro_smarts, popularity=len(rxns))
        rec['_id'] = rec_count.next()
        jsonschema.validate(rec, SCHEMA)
        records.append(rec)

    client = MongoClient()
    db = client[args.database]
    coll = db[args.collection]
    for rec in records:
        coll.insert(rec)


def get_transforms(filename):
    with open(filename) as f:
        data = pickle.load(f)
    return data


def reverse_smarts(smarts):
    reactant_templates, product_templates = smarts.split('>>')
    return '>>'.join([product_templates, reactant_templates])


def make_record(smarts, popularity=0):
    rec = {
        'reaction_smarts': smarts,
        'popularity': popularity,
    }
    return rec


if __name__ == '__main__':
    main()