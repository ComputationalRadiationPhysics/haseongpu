#!/bin/sh


if [ "$1" = "" ] ; then
	echo "Usage: $0 <folder/filename.mat>"
	exit 0
fi

INPUT="$1"
if [ "$(echo "$INPUT" | grep -c '.*.mat')" = "0" ] ; then
	echo "$INPUT" is no valid *.mat file!
	exit 1
fi

if [ $(echo "$INPUT" | grep -c -e ".*/") -ne 1 ] ; then
	echo 'You did not specify a folder!'
	exit 1
fi

FOLDER="$(echo $INPUT | grep -o -e '.*/')"
FILE="$(echo $INPUT | grep -o -e '/.*.mat' | cut -c 2-)"


rm -f ${FOLDER}*.txt
rm -f ${FOLDER}raw_input.zip
octave --silent --eval "parse_mat('$FOLDER','$FILE')" 1>/dev/null

cd $FOLDER
zip raw_input.zip *.txt 1>/dev/null
rm -f *.txt
