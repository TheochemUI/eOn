Serve mode implements `Potential.getCapabilities` on the single callback server and the pooled server.
The reply names protocol family `rgpot.potentials`, protocol major 1, schema id `0xbd1f89fa17369103`, bridge ABI major 1 and minor 0, layout 1, DLPack 1.0, operations energy and forces, and `bridgeFeatures` 0.
`buildVersion` is the project version. `buildRevision` is `git rev-parse --short=12`, or empty when that command fails.
A client that calls `getCapabilities` before `calculate` refuses a server that does not implement the method.
