On Windows the metatomic engine test loads torch's `uv.dll` from the torch install. The pixi prefix has another `uv.dll`, and that copy does not provide the entry point torch was linked against.
