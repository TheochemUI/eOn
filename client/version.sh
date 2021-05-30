# Todo, handle this more gracefully
VERSION=$(svnversion 2> /dev/null)
if [ -z "$VERSION" ]; then
    VERSION=$(git rev-parse --short HEAD 2> /dev/null)
fi


if [ -z "$VERSION" ]; then
    VERSION="unknown"
fi

echo static const char VERSION[] = \"$VERSION\"\; 
echo static const char BUILD_DATE[] = \"`date`\"\; 
