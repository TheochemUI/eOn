
#include "Bundling.h"

#include <cstdio>
#include <dirent.h>
#include <errno.h>
#include <cstring>
#include <cstdlib>
#include <unistd.h>

int copyfile(char *, char *);

int getBundleSize(void)
{
    DIR *dir;
    struct dirent *dp;
    int num_bundle=-1;

    dir = opendir(".");

    //Find the highest numbered bundle file
    //and return that number.
    while ((dp=readdir(dir))) {
        if (dp->d_name[0] == '.') {
            continue;
        }

        //If "config" is not in the filename
        //then skip.
        if (strstr(dp->d_name, "config")==NULL &&
            strstr(dp->d_name, "ini")==NULL) {
            continue;
        }

        //Find the last underscore
        char *ch = strrchr(dp->d_name, '_');
        if (ch == NULL) {
            continue;
        }else{
            ch += 1;
        }
        //Find the last period
        char *cch = strrchr(dp->d_name, '.');
        if (cch == NULL) continue;
        *cch = '\0';
        if (isdigit(*ch)) {
            int i=atoi(ch)+1;
            if (i>num_bundle) {
                num_bundle = i;
            }
        }
    }

    closedir(dir);

    return num_bundle;
}

int strchrcount(char *haystack, char needle)
{
    int count=0;
    for (char *ch=haystack;*ch!='\0';ch++) {
        if (*ch == needle) {
            count++;
        }
    }
    return count;
}

std::vector<std::string> unbundle(int number) {
    std::vector<std::string> filenames;
    DIR *dir;
    struct dirent *dp;

    dir = opendir(".");

    while ((dp=readdir(dir))) {
        int bundleNumber;
        std::string originalFilename(dp->d_name);

        if (dp->d_name[0] == '.') {
            continue;
        }

        //If "passed" is not in the filename
        //then skip.
        //if (strstr(dp->d_name, "passed")==NULL) {
        //    continue;
        //}

        int numUnderscores = strchrcount(dp->d_name, '_');
        if (numUnderscores < 1) {
            continue;
        }

        //Find the last underscore
        char *ch = strrchr(dp->d_name, '_');
        if (ch == NULL) {
            continue;
        }else{
            ch += 1;
        }
        //Find the last period
        char *cch = strrchr(dp->d_name, '.');
        if (cch == NULL) continue;
        *cch = '\0';
        if (isdigit(*ch)) {
            bundleNumber=atoi(ch);
            if (bundleNumber != number) {
                continue; 
            }
        }

        *(ch-1) = '\0';

        int stringSize = strlen(dp->d_name) + 10;
        char *newFilename = new char[stringSize];
        snprintf(newFilename, stringSize, "%s.%s", dp->d_name, cch+1);

        int err;
        err = copyfile((char *)originalFilename.c_str(), newFilename);
        if (err) {
            fprintf(stderr, "error: unbundle: problem copying %s to %s\n",
                    originalFilename.c_str(), newFilename);
        }
        filenames.push_back(std::string(newFilename));
        delete[] newFilename;
    }

    closedir(dir);
    return filenames;
}

//Why couldn't there be a crossplatform standard for file copy?
int copyfile(char *filenameSrc, char *filenameDest)
{
    //10kB buffer
    int buffSize = 1024*10;
    char *buff = new char[buffSize];
    FILE *fSrc = fopen(filenameSrc, "rb"); 
    if (fSrc == NULL) {
        fprintf(stderr, "error: copyfile: problem opening src file\n");
        delete [] buff;
        return 1;
    }
    FILE *fDest = fopen(filenameDest, "wb");
    if (fDest == NULL) {
        fprintf(stderr, "error: copyfile: problem opening dest file\n");
        delete [] buff;
        return 1;
    }

    //XXX: Error checking
    while (!feof(fSrc)) {
        int count = fread(buff, sizeof(char), buffSize, fSrc);
        fwrite(buff, sizeof(char), count, fDest);
    }
    
    delete [] buff;
    fclose(fSrc);
    fclose(fDest);

    return 0;
}

void deleteUnbundledFiles(std::vector<std::string> unbundledFilenames)
{
    for (unsigned int i=0;i<unbundledFilenames.size();i++) {
        unlink(unbundledFilenames[i].c_str());
    }
}

void bundle(int number, std::vector<std::string> filenames,
            std::vector<std::string> *bundledFilenames) {
    for (unsigned int i=0;i<filenames.size();i++) {

        std::string filename = filenames[i];
        //Assumes that bundle numbers wont exceed 100,000,000
        int stringSize = filename.length()+20;
        char *newFilename;
        newFilename = new char[stringSize];
        strncpy(newFilename, filename.c_str(), stringSize);

        char *ch;
        char *fileEnding;
        if ((ch=strrchr(newFilename, '.')) != NULL) {
            *ch = '\0';
            fileEnding = new char[filename.length()];
            strncpy(fileEnding, ch+1, filename.length());
            char *buff = new char[stringSize];
            snprintf(buff, stringSize, "%s_%i.%s", 
                     newFilename, number, fileEnding);
            strncpy(newFilename, buff, stringSize);
            delete [] fileEnding;
            delete [] buff;
        }else{
            snprintf(newFilename, stringSize, "%s_%i",
                     newFilename, number);
        }

        int err;
        if ((err=rename(filename.c_str(), newFilename))) {
            fprintf(stderr, "error: bundle: cannot rename %s to %s: %s\n",
                    filename.c_str(), newFilename, strerror(errno));
        }
        bundledFilenames->push_back(std::string(newFilename));

        delete [] newFilename;
    }
}
