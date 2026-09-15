IRACompare.cpp is always compiled so NEB can call match_endpoints when
IRA is off (the methods stub). `matchArrays` plus `pyeonclient.ira_match`
let `eon.atoms.rot_match` use IRA when the wheel has it, else Kabsch.
