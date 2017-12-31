#!/bin/bash

# Find job IDs
foo=`qstat -u $USER | grep -oE '^[0-9]+'`