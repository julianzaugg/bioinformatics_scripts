#!/usr/bin/env bash

# Takes the output from bbtools readlength.sh and calculates the gigabases

cat $1 | grep "^[^#;]" | awk '{print $0, $1*$2/1000000000}' | awk '{c+=$10}END{print c}'
