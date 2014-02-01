#!/bin/bash

FILES="*.dfx"

for i in $FILES
do
    dsr -i $i
done




