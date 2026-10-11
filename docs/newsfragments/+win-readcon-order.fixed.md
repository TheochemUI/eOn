A Windows wheel builds readcon-core before linking it. The wrap fallback had dropped the order stamp, so link.exe opened readcon_core.lib before cargo wrote the file.
