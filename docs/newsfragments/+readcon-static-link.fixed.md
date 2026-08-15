The readcon wrap overlay builds only the static archive as a single
custom_target output, so Meson waits for libreadcon_core.a before
linking eonclib and shared objects do not record NEEDED
libreadcon_core.so.