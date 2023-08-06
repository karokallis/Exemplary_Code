#!/bin/bash

cd $1

for i in *
do
 echo "$i"
 output=$(echo "$i" | sed 's/..$//')

 mv "$i" "$output"
done
