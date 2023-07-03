#!/usr/bin/env python

from optparse import OptionParser
import sys

def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def relative_error(x,y):
    if x==y:
        return 0.0
    elif x != 0.0 and y != 0.0:
        return abs(x-y)/min(abs(x),abs(y))
    elif y == 0.0 and x != 0.0:
        return 1.0
    elif x == 0.0 and y != 0.0:
        return 1.0
    else:
        return 0.0

def ndiff(filename1, filename2, rel_tol=0.0, quiet=True):
    f1 = open(filename1).readlines()
    f2 = open(filename2).readlines()

    same = True
    max_rel_err = 0.0
    reason = ""

    if len(f1) != len(f2):
        reason = "files are of different length"
        return False, max_rel_err, reason

    for i in range(len(f1)):
        r1 = f1[i].split()
        r2 = f2[i].split()

        if len(r1) != len(r2):
            reason = "line %i has a different number of records" % (i+1)
            return False, max_rel_err, reason

        for j in range(len(r1)):
            r1_isnum=False
            r2_isnum=False
            if is_number(r1[j]):
                r1[j] = float(r1[j])
                r1_isnum = True
            if is_number(r2[j]):
                r2[j] = float(r2[j])
                r2_isnum = True

            if r1_isnum != r2_isnum:
                reason = "record %i on line %i are not of matching types" % (j+1,i+1)
                return False, max_rel_err, reason

            if r1_isnum == False and r2_isnum == False:
                if r1[j] != r2[j]:
                    reason = "record %i on line %i are not matching strings" % (j+1,i+1)
                    return False, max_rel_err, reason
            else:
                rel_err = relative_error(r1[j],r2[j])

                if rel_err > max_rel_err:
                    max_rel_err = rel_err

                if rel_err > rel_tol:
                    if not quiet:
                        print("relative error of %.3e of record %i on line %i" % (rel_err,j+1,i+1))
                    same = False
                    reason = "maximum relative error exceeds tolerance %.3e" % rel_tol

    return same, max_rel_err, reason

def main():
    usage = "usage: %prog [options] filename1 filename2"
    parser = OptionParser(usage=usage)

    parser.add_option("-r", "--relative-error", dest="relerr", default=0.0,
                      help="the maximum relative error")
    parser.add_option("-q", "--quiet",
                      action="store_true", dest="quiet", default=False,
                      help="don't print a line per numerical difference found")
    parser.add_option("-s", "--silent", action="store_true", dest="silent",
                      default=False, help="don't print anything"
                      " this is useful if you just want to use the exit status")

    (options, args) = parser.parse_args()

    if len(args) != 2:
        print("too few arguments")
        parser.print_usage()
        sys.exit(2)

    if options.silent:
        options.quiet = True

    same = ndiff(args[0], args[1], options.relerr, options.quiet)

    if same:
        sys.exit(0)
    else:
        sys.exit(1)

if __name__ == '__main__':
    main()
