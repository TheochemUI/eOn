import subprocess
from os.path import abspath, dirname


def version():
    path = dirname(abspath( __file__ ))
    try:
        output = subprocess.check_output(["svnversion",str(path)], stderr=subprocess.STDOUT)
    except subprocess.CalledProcessError as grepexc:                                                                                                   
        return 'unknown'
    return 'svn revision %s' % output.decode('ascii')


if __name__ == '__main__':
    print(version())
