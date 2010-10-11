#include <errno.h>
#include <string.h>
#include <stdio.h>
#include <fcntl.h>
#include <dirent.h>
#include <string.h>
#include <sys/stat.h>

#include <archive.h>
#include <archive_entry.h>

#include "Compression.h"

//const long potentialNewPotential = 1;

//Probably should be a multiple of 512 as tar files have a
//block size of 512 bytes.
#define BUFFER_SIZE  10240

// This function will create a gzipped tar file (outname)
// which contains the files matched by the function pattern_match.
// Pattern_match should return 0 on no match and 1 on match.
int create_archive(char *outname, char *path, int (*pattern_match)(char *))
{
    struct archive *a;
    struct archive_entry *entry;
    struct stat st;
    char buff[BUFFER_SIZE];
    int len;
    int fd;
    DIR *dp;
    struct dirent *ep;

    a = archive_write_new();
    archive_write_set_compression_gzip(a);
    archive_write_set_format_pax_restricted(a);
    archive_write_open_filename(a, outname);

    dp = opendir(path);
    if (dp != NULL) {
        while ((ep = readdir(dp))) {
            /* This had to be removed for win32 compat :(
            if (ep->d_type != DT_REG) {
                continue;
            }
            */
            if (pattern_match(ep->d_name)==0) {
                continue;
            }
            stat(ep->d_name, &st);
            entry = archive_entry_new();
            archive_entry_set_pathname(entry, ep->d_name);
            archive_entry_set_size(entry, st.st_size);
            archive_entry_set_filetype(entry, AE_IFREG);
            archive_entry_set_perm(entry, 0644);
            archive_write_header(a, entry);
            fd = open(ep->d_name, O_RDONLY);
            len = read(fd, buff, sizeof(buff));
            while (len > 0) {
                archive_write_data(a, buff, len);
                len = read(fd, buff, sizeof(buff));
            }
            close(fd);
            archive_entry_free(entry);
        }
        closedir(dp);
    }else{
        printf("error opening directory\n");
    }

    archive_write_close(a);
    archive_write_finish(a);

    return 0;
}

int extract_archive(char *filename)
{
    struct archive *a;
    struct archive_entry *entry;
    int r;
    size_t size;
    char buff[BUFFER_SIZE];
    FILE *fd;

    a = archive_read_new();
    archive_read_support_compression_gzip(a);
    archive_read_support_format_tar(a);

    r = archive_read_open_filename(a, filename, 10240);

    if (r != ARCHIVE_OK) {
        return 1;
    }

    for (;;) {
        if ((r = archive_read_next_header(a, &entry))) {
            if (r != ARCHIVE_OK) {
                if (r == ARCHIVE_EOF) {
                    return 0;
                }else{
                    return 1;
                }
            }
        }

        fd = fopen(archive_entry_pathname(entry),"wb");
        if (fd == NULL) {
            fprintf(stderr, "problem extracting archive: %s: %s\n", filename,
                    strerror(errno));
            return 1;
        }
        for (;;) {
            size = archive_read_data(a, buff, BUFFER_SIZE);
            if (size < 0) {
                return 1;
            }

            if (size == 0) {
                break;
            }

            fwrite(buff, 1, size, fd);
        }
        fclose(fd);
    }

    r = archive_read_finish(a);
    if (r != ARCHIVE_OK) {
        return 4;
    }
    return 0;
}
