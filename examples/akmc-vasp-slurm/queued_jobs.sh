#!/bin/sh
sqs | awk '$1~/^[0-9]/ {print $1}'
