#!/bin/bash


LATEX_MAIN_DOCUMENT="talk"


red='\e[0;31m'
green='\e[0;32m'
yellow='\e[0;33m'
NC='\e[0m' #no color

IGNORE1="Package: infwarerr 2010/04/08 v1.3 Providing info/warning/error messages (HO)"
IGNORE2="Class scrreprt Warning: Usage of package \`fancyhdr' together"



function wait_fixed(){
  echo "Press key to try again!"
  read
}

if [ "$1" != "" ] ; then
  LATEX_MAIN_DOCUMENT="$1"
fi

while [ 1 ] ; do

  make
  if [ $? -ne 0 ] ; then
    wait_fixed
    continue
  fi

  clear
  echo "Last compilation done at `date`"
  N_WARNINGS=$(grep -i warning "$LATEX_MAIN_DOCUMENT".log |  grep -c -v -e "$IGNORE1" -e "$IGNORE2" )
  if [ "$N_WARNINGS" -gt 0 ] ; then
    echo -e "\n${yellow}Warnings:\n"
    grep -i warning "$LATEX_MAIN_DOCUMENT".log |  grep -v -e "$IGNORE1" -e "$IGNORE2"
    echo -e "${NC}\n\n\n\n"
  fi

  N_ERRORS=$(grep -i error "$LATEX_MAIN_DOCUMENT".log |  grep -c -v -e "$IGNORE1" -e "$IGNORE2" )
  if [ "$N_ERRORS" -gt 0 ] ; then
    echo -e "\n${red}Errors:\n"
    grep -i error "$LATEX_MAIN_DOCUMENT".log |  grep -v -e "$IGNORE1" -e "$IGNORE2"
    echo -e "${NC}\n\n\n\n"
  fi


  if [ "$N_WARNINGS" -eq 0 ] && [ "$N_ERRORS" -eq 0 ] ; then
    echo -e "${green}"
    echo '       /(|'
    echo '      (  :'
    echo '     __\  \  _____'
    echo '   (____)  `|'
    echo '  (____)|   |'
    echo '   (____).__|'
    echo '    (___)__.|_____'
    echo -e "${NC}"
  fi 

  

	notify-send -t 3000 "pdflatex" "compilation done"

	inotifywait -q -e modify *.tex slides/*.tex
done


