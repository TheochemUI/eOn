#!/usr/bin/env bash

VERSION=$(svnversion 2> /dev/null)

if [ -z "$VERSION" ]; then
    VERSION="unknown"
fi

echo static const char VERSION[] = \"$VERSION\"\; 
echo static const char BUILD_DATE[] = \"`date`\"\; 
