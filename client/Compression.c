#include <stdio.h>
#include <fcntl.h>
#include <sys/stat.h>

#include <archive.h>
#include <archive_entry.h>

#include "compression.h"

//Probably should be a multiple of 512 as tar files have a
//block size of 512 bytes.
#define BUFFER_SIZE  10240

/*
int main(int argc, char *argv[])
{
    int r;
    if ((r = extract_archive("input.dat")) != 0) {
        printf("error opening archive: %i\n", r);
    }
    const char *filenames[] = {"compression.c", "compression.h"};
    if ((r = create_archive("output.dat", filenames)) != 0) {
        printf("error opening archive: %i\n", r);
    }
}
*/

// This function will create a gzipped tar file (outname)
// which contains the files in **filename
int create_archive(const char *outname, const char **filename)
{
    struct archive *a;
    struct archive_entry *entry;
    struct stat st;
    char buff[BUFFER_SIZE];
    int len;
    int fd;

    a = archive_write_new();
    archive_write_set_compression_gzip(a);
    archive_write_set_format_pax_restricted(a);
    archive_write_open_filename(a, outname);

    while (*filename) {
        stat(*filename, &st);
        entry = archive_entry_new();
        archive_entry_set_pathname(entry, *filename);
        archive_entry_set_size(entry, st.st_size);
        archive_entry_set_filetype(entry, AE_IFREG);
        archive_entry_set_perm(entry, 0644);
        archive_write_header(a, entry);
        fd = open(*filename, O_RDONLY);
        len = read(fd, buff, sizeof(buff));
        while (len > 0) {
            archive_write_data(a, buff, len);
            len = read(fd, buff, sizeof(buff));
        }
        close(fd);
        archive_entry_free(entry);
        filename++;
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
    archive_read_support_compression_all(a);
    archive_read_support_format_all(a);

    r = archive_read_open_filename(a, filename, BUFFER_SIZE);

    if (r != ARCHIVE_OK) {
        return 1;
    }

    for (;;) {
        if (r = archive_read_next_header(a, &entry)) {
            if (r != ARCHIVE_OK) {
                if (r == ARCHIVE_EOF) {
                    return 0;
                }else{
                    return 1;
                }
            }
        }

        fd = fopen(archive_entry_pathname(entry),"wb");
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
