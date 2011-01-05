#!/bin/bash
# Compiles Xwalk from source code.
# Requires SUN java 1.6.
if [ ! -d bin ]
then
    mkdir bin
fi

# compile source files and add binaries to bin/ directory
javac -sourcepath src/ -d bin/ src/Xwalk.java

# add parameter text files to bin/ directory
mkdir -p bin/mm/parameter
cp src/mm/parameter/* bin/mm/parameter
