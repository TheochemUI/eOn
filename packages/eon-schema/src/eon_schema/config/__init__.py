"""L1 full job-config models (pydantic).

These are the eOn configuration models used for docs validation and as the
shared schema for **eon-akmc** and any consumer that needs the full job
graph. They live in this package so L0 (Cap'n Proto) and L1 share one home.

Import paths (equivalent)::

    from eon_schema.config import MainConfig, Config
    from eon.schema import MainConfig, Config  # re-export in eon-akmc

L2 request models (``DimerSpec``, ``NebSpec``) stay under ``eon_schema.api``.
"""

from __future__ import annotations

from eon_schema.config.models import *  # noqa: F403
from eon_schema.config.models import (  # noqa: F401
    MainConfig,
    StructureComparisonConfig,
    AKMCConfig,
    BasinHoppingConfig,
    PathsConfig,
    CommunicatorConfig,
    ProcessSearchConfig,
    PrefactorConfig,
    PotentialConfig,
    XTBPot,
    ZBLPot,
    SocketNWChemPot,
    RgpotPot,
    Metatomic,
    ASE_NWCHEM,
    ASE_ORCA,
    AMSConfig,
    AMSIOConfig,
    AMSEnvConfig,
    SaddleSearchConfig,
    KDBConfig,
    RecyclingConfig,
    CoarseGrainingConfig,
    OptimizerConfig,
    QuickMinConfig,
    FIREConfig,
    LBFGSConfig,
    CGConfig,
    SDConfig,
    XtsciConfig,
    RefineConfig,
    ServeConfig,
    DebugConfig,
    DimerConfig,
    NudgedElasticBandConfig,
    LanczosConfig,
    DavidsonConfig,
    ARTnConfig,
    IRAConfig,
    BGSDConfig,
    HessianConfig,
    DynamicsConfig,
    ParallelReplicaConfig,
    HyperdynamicsConfig,
    GPRDimerConfig,
    DistributedReplicaConfig,
    Config,
)

__all__ = [
    "MainConfig",
    "StructureComparisonConfig",
    "AKMCConfig",
    "BasinHoppingConfig",
    "PathsConfig",
    "CommunicatorConfig",
    "ProcessSearchConfig",
    "PrefactorConfig",
    "PotentialConfig",
    "XTBPot",
    "ZBLPot",
    "SocketNWChemPot",
    "RgpotPot",
    "Metatomic",
    "ASE_NWCHEM",
    "ASE_ORCA",
    "AMSConfig",
    "AMSIOConfig",
    "AMSEnvConfig",
    "SaddleSearchConfig",
    "KDBConfig",
    "RecyclingConfig",
    "CoarseGrainingConfig",
    "OptimizerConfig",
    "QuickMinConfig",
    "FIREConfig",
    "LBFGSConfig",
    "CGConfig",
    "SDConfig",
    "XtsciConfig",
    "RefineConfig",
    "ServeConfig",
    "DebugConfig",
    "DimerConfig",
    "NudgedElasticBandConfig",
    "LanczosConfig",
    "DavidsonConfig",
    "ARTnConfig",
    "IRAConfig",
    "BGSDConfig",
    "HessianConfig",
    "DynamicsConfig",
    "ParallelReplicaConfig",
    "HyperdynamicsConfig",
    "GPRDimerConfig",
    "DistributedReplicaConfig",
    "Config",
]

# INI helpers (config write / hydrate without eon-akmc)
from eon_schema.config.ini import (  # noqa: E402,F401
    INI_FIELD_ALIASES,
    MODEL_INI_SECTION,
    allowed_keys,
    defaults_from_catalog,
    format_ini_value,
    hydrate_ini,
    model_to_ini_section,
    models_to_ini,
    read_ini,
    unknown_ini_keys,
    write_ini,
    write_models_ini,
)

__all__ += [
    "INI_FIELD_ALIASES",
    "MODEL_INI_SECTION",
    "allowed_keys",
    "defaults_from_catalog",
    "format_ini_value",
    "hydrate_ini",
    "model_to_ini_section",
    "models_to_ini",
    "read_ini",
    "unknown_ini_keys",
    "write_ini",
    "write_models_ini",
]
