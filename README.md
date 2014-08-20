# Retronet

A package allowing to generate a retrosynthetic chemical network around a
specified target.


## Requirements

To use `retronet` package, the following third-party packages must be installed
on your system:

*	[RDKit](http://www.rdkit.org): open source toolkit for chemoinformatics;
*	[NetworkX](https://networkx.github.io): high-productivity software for
	complex networks.

Also, retrosynthetic transforms, encoded in
[Daylight's](http://www.daylight.com/) SMARTS notation, must be stored as a
[mongodb](http://www.mongodb.org/) collection in a location accessible from
your box.

The detailed description of a transform database entry, using
[JSON-schema](http://json-schema.org/), can be found in `transform.json`.


## Usage

At the minimum, running the program requires providing a text file containing
chemical targets encoded in [Daylight's](http://www.daylight.com/) SMILES
notation (eg. 'CCO' for ethyl alcohol), i.e.,

	python2 netgen.py targets.smi

where `targets.smi` contains SMILES of desired targets, each in a separate
line.  For full list of supported options, type

    python2 netgen.py --help
