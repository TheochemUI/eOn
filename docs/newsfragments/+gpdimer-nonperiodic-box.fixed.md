`AtomicGPDimer` passes a zero cell into the GP force box when
`Matter::getPeriodic` is false. XTB no longer treats a stored .con box
as periodic boundaries on those force calls.
