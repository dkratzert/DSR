#!/bin/bash

# initialize DSR_DIR to default directory if not set properly
if [ ! $DSR_DIR ]; then 
    export DSR_DIR=/opt/DSR
fi

PATH=$PATH:$DSR_DIR
export PATH
export DSR_DIR