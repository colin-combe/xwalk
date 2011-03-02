if [ ! -d bin ]
    then  
    mkdir bin
fi
if [ ! -d bin/mm/parameter ] 
    then
    mkdir -p bin/mm/parameter
fi
if [ ! -d bin/xwalk/parameter ] 
    then
    mkdir -p bin/xwalk/parameter
fi
cp -f src/mm/parameter/*.txt bin/mm/parameter/
cp -f src/xwalk/parameter/*.txt bin/xwalk/parameter/
javac -sourcepath src/ -d bin/ src/Xwalk.java
