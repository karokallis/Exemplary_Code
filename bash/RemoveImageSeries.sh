#!/bin/bash

#InputPath = $1
cd $1

for f in */; do
    if [ -d "$f" ]; then
        #echo "${f}IMG00001.dcm"
        dcm=$(find ${f} -type f | head -n 1)
        #echo ${dcm}
        sd=$(/usr/pubsw/packages/dcmtk/3.6.0/bin/dcmdump ${dcm} | grep -i "(0008,103e)");
        if [[ ${sd} == *"RSI"* || ${sd} == *"T2"*]]; then
            echo ${sd}
        fi
        
        #sd=$(dcmdump ${dcm}
        #echo ${sd}
        #for dcm in "${f}"
        #do
            #echo ${dcm}
            #sd=$(dcmdump ${dcm} | grep 'SeriesDescription');
        #sd = $(dcmdump "$f$VAR2" | grep 'SeriesDescription'
        #done
    fi
done