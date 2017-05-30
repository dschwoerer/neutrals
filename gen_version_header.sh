#!/bin/bash

echo "Generating header for git hash"
#$1 $2
#GIT_HEADER="$1/$2"
#if [ -z "$2" ]; then
GIT_HEADER="./git_version.hxx"
#fi

GIT_VERSION=`git rev-parse HEAD`
#GIT_VERSION="`git -C \"$1\" describe`"
if grep --quiet $GIT_VERSION $GIT_HEADER; then
    echo "No need to generate new $GIT_HEADER - git hash is unchanged"
    exit 0;
fi

echo "git version is:" $GIT_VERSION

echo "#pragma once" > $GIT_HEADER
echo "" >> $GIT_HEADER
echo "const char* NEUTRALS_GIT_SHA1 = \"$GIT_VERSION\";" >> $GIT_HEADER
echo >> $GIT_HEADER
echo "file is generated into" $GIT_HEADER
