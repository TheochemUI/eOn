The in-process path keeps the typed job result and does not build a results.dat buffer. Cluster adapters still format that text, and parse_results stays the reader for file workers.
