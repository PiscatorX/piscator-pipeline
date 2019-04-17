#!/bin/bash

grep  ">" "${1}"  | tr -d ">" | cut -d " " -f 1,2  --output-delimiter "	"
