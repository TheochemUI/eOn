Windows pyeonclient wheels build the client libraries as static archives, matching the Windows CI build. A plugin DLL with no exports never wrote the import library MSVC then failed to open.
