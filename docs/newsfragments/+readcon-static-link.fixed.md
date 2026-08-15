The readcon wrap overlay still cargo-builds both crate types. A generated
.c that depends on that custom_target is the ninja order edge, so eonclib
does not link before the outputs exist and does not record NEEDED
libreadcon_core.so or pass the DLL to MSVC.