#!/bin/bash
CMD=$1
FILE=$2
log_postfix=$3
echo "starting" > "$FILE"_"$log_postfix"
eval "$CMD"
echo "done" > "$FILE"_"$log_postfix"
