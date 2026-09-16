from collections import Counter
from pathlib import Path

import eon

i = 0
while (Path("states") / str(i)).is_dir():
	p = eon.fileio.loadcon(str(Path("states") / str(i) / "reactant.con"))
	cna = eon.atoms.cnat(p,4.0)
	count = Counter()
	for c in cna:
		count[int(c)] += 1
	tally = []
	for k, v in count.items():
		tally.append((k, v))
	tally.sort(key = lambda k: k[0])
	for t in tally:
		print("%3d" % t[1], end=" ")
	print()

	for j in range(len(p)):
		p.names[j] = eon.atoms.elements[int(cna[j]) + 1]['symbol']
	eon.fileio.savecon('movie.con', p, w='a')

	i += 1
