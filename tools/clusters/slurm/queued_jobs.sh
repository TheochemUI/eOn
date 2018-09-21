#!/bin/sh
squeue -u graeme | awk '$1~/^[0-9]/ {print $1}'
