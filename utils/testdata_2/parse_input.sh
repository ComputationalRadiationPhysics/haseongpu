#!/bin/sh

rm *.txt
rm raw_input.zip
octave parse_mat.m 1>/dev/null 
zip raw_input.zip *.txt 1>/dev/null
rm *.txt
