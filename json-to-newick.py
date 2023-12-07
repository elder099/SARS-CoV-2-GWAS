from augur.utils import json_to_tree
import Bio.Phylo
import json

# Opening JSON file
f = open(os.path.expanduser('~/CDPH/GWAS/GWASv6.0/output/nextclade.auspice.json'))

# returns JSON object as
# a dictionary
tree_json = json.load(f)

# Convert JSON to BioPython tree instance.
tree = json_to_tree(tree_json)

# Write tree to disk as a Newick file.
Bio.Phylo.write(tree, os.path.expanduser('~/CDPH/GWAS/GWASv6.0/Newicktree.nwk'), "newick")
