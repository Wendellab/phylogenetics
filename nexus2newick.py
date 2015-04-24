#convert nexus tree format to newick

import sys
import os
import dendropy

for n in sys.argv[1:]:
    basename = os.path.splitext(n)
    outfile = basename[0] + ".nwk"
    nexusfile = dendropy.TreeList.get_from_path(n, "nexus")
    nexusfile.write_to_path('temp.nwk', "newick")
    os.rename('temp.nwk', outfile)
