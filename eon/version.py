import subprocess
from os.path import abspath, dirname


# TODO(rg): Grab a dynamic version..
def version():
    path = dirname(abspath( __file__ ))


if __name__ == '__main__':
    print(version())
