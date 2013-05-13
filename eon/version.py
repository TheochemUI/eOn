import commands
from os.path import abspath, dirname


def version():
    path = dirname(abspath( __file__ ))
    status, output = commands.getstatusoutput('svnversion %s' % path) 
    if status != 0:
        return 'unknown'
    else:
        return 'svn revision %s' % output

if __name__ == '__main__':
    print version()
