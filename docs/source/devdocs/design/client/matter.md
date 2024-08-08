# Matter design

The design of Matter is meant to be orthogonal to the reading and writing of
different outputs. This is conceptually close to the design of `Potential`.

We need to set some boundaries about what Matter is and what it is not meant to
be. In practice, it is desirable to be able to construct empty matter objects,
especially because having to correctly initialize `Matter` before use can be
restrictive, i.e. lazy loading of data is desirable.

## Populating

```{versionchanged} 3.x
- Matter objects do not hold headers, these are stored in ConFileParser
```


```{versionadded} 3.x
```

- `Matter` is now populated by `ConFileParser`, which handles validation and works for both `con` and `convel` files .

## Caching

```{versionadded} 3.x
```

Earlier versions of `eOn` used the concept of a tainted bit to provide a way to
short-circuit recalculations at the same state. Essentially, the
`recomputePotential` boolean was set for certain operations and controlled when
/ how the potential energy and forces were computed. This was essentially an LRU
cache of size one, which notably prevented const-correctness of methods, and
prevented re-use of previous calculations.

The new design requires `Matter` to be constructed with a `cachelot` instance,
with the parameters set at runtime but offers the ability to re-use older
calculations and also for different `Matter` objects to share the same cache. 

```{code-block} cpp
const auto config = toml::table{{"Potential", toml::table{{"potential", "LJ"}}}};

std::string confile = "pos.con";
eonc::mat::ConFileParser cfp;
auto pot1 = makePotential(config);
auto pot2 = makePotential(config);
auto pot3 = makePotential(config);

auto CACHELOT_EONCTEST = cachelot::cache::Cache::Create(eonc::cache_memory,
eonc::page_size, eonc::hash_initial, true);

Matter mat1(pot1, &CACHELOT_EONCTEST);
cfp.parse(mat1, confile);

Matter mat2(pot2, &CACHELOT_EONCTEST);
cfp.parse(mat2, confile);

Matter mat3(pot3, &CACHELOT_EONCTEST);
cfp.parse(mat3, confile);

mat1.getPotentialEnergy();
REQUIRE(pot1->getTotalForceCalls() == 1);

mat2.getPotentialEnergy();
REQUIRE(pot2->getTotalForceCalls() == 1);

mat3.getPotentialEnergy();
REQUIRE(pot3->getTotalForceCalls() == 1);
REQUIRE(pot1->getInstances() == 3);
```

```{todo}
```

Eventually:

- Run in an async thread (TODO(rg))
- Be used for structure comparison (TODO(rg))
- Allow on-the-fly training of potentials for positions close to existing ones (TODO(rg))
