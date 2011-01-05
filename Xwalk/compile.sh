#!/bin/bash
# requires SUN java 1.6
if [ ! -d bin ]
then
    mkdir bin
fi

SRC=$(find | grep ".java$")
javac -d bin $SRC
