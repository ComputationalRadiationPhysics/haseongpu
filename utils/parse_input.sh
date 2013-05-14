#!/bin/sh


if [ "$1" = "" ] ; then
	echo "Usage: $0 <folder/filename.mat>"
	exit 0
fi

INPUT="$1"
if [ $(echo "$INPUT" | grep -c *.mat) -ne 1 ] ; then
	echo '"$INPUT" is no valid *.mat file!'
	exit 1
fi

if [ $(echo "$INPUT" | grep -c -e ".*/") -ne 1 ] ; then
	echo 'You did not specify a folder!'
	exit 1
fi

FOLDER="$(echo $INPUT | grep -e '.*/')"
FILE="$(echo $INPUT | grep -e '/.*.mat' | cut -c 2-)"

cd $FOLDER

rm *.txt
rm raw_input.zip
octave -qf --eval parse_mat("$FILE")
zip raw_input.zip *.txt 1>/dev/null
rm *.txt
