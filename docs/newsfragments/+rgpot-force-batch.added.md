`RgpotAdapter` reports `caps().batched` and forwards a band of images to
`rgpot::forceBatchImpl`, so a kernel that evaluates several systems in
one call sees the whole band. The wrap tracks the rgpot branch that
adds `ForceBatch` until a tag carries it.
