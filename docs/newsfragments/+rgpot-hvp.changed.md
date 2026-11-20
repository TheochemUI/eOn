Pin the rgpot wrap to `758dbc02`, the merge of OmniPotentRPC/rgpot#292 into
`develop`, which also contains the analytic Hessian-vector products merged in
#258. RgpotAdapter forwards those products as eOn mass-weighted
Hessian-vector products. Existing InstantonJob calculations use the exact
product when available instead of a potentially step-size-biased
finite-difference Hessian. Remove eOn's vendored rgpot Cap'n Proto include
patch because upstream now propagates the required RPC header dependency.
