The readcon wrap overlay attaches only the static custom_target output
as a source, so Meson still builds the archive before link and shared
objects do not record NEEDED libreadcon_core.so.