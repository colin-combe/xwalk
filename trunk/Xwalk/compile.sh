#!/bin/bash
# requires SUN java 1.6
if [ ! -d bin ]
then
    mkdir bin
fi
javac -sourcepath src/ -d bin/ src/Xwalk.java
