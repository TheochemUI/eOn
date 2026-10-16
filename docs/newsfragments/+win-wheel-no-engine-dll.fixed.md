The Windows pyeonclient wheel no longer installs eon_relax_engine.dll. meson-python will not bundle a Windows DLL, and the extension does not load that module.
