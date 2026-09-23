Langevin dynamics no longer moves a coordinate frozen on only one axis.
`Dynamics::langevinVerlet` applies noise and the position step per free axis, not only when the whole atom is fixed.
