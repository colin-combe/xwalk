#!/bin/bash
if [ ! -d bin ]
then
    mkdir bin
fi
javac -sourcepath src/ -d bin/ src/Xwalk.java
