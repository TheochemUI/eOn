---
myst:
  html_meta:
    "description": "The client runs optimal hyperplane transition state theory."
    "keywords": "eOn hyperplane, free energy, pmf scan"
---

# Optimal hyperplane transition state theory

The method places a dividing plane between the reactant and the product.
The client token is `oh_tst`. The rendered job list omits that token.
A schema check rejects the name. The client accepts the name in
`config.ini`.

The section is `[OH_TST]`. `reactant_filename` defaults to `pos.con`.
`product_filename` defaults to `product.con`. `equil_steps` defaults to
200. `sample_steps` defaults to 800. Both counts are the constrained
sampling length on each plane. `max_planes` defaults to 200.
`thermostat` is `andersen` (the default) or `gle`. `symmetry_products`
names other product files, separated by commas. `pmf_scan` keeps the
plane normal on the reactant-to-product line. It places `scan_planes`
(default 40) uniformly along that line.

```{code-block} ini
[Main]
job = oh_tst
temperature = 300

[OH_TST]
reactant_filename = pos.con
product_filename = product.con
equil_steps = 200
sample_steps = 800
max_planes = 200
```
