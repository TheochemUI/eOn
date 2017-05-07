#!/bin/sh
qstat | awk '$1~/^[0-9]/ {print $1}'
