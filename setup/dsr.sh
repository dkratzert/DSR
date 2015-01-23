#!/bin/bash

# initialize DSR_DIR to default directory if not set properly
if [ ! $DSR_DIR ]; then 
    export DSR_DIR=/opt/DSR
fi

if [ ! $DSR_DB_DIR ]; then 
    export DSR_DB_DIR=/opt/DSR
fi

PATH=$PATH:$DSR_DIR
DSR_DIR=$DSR_DB_DIR
export PATH
export DSR_DB_DIR
export DSR_DIR