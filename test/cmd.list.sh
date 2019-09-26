#!/bin/bash
#$ -S /bin/bash -tc 1
if [ $SGE_TASK_ID -eq 1 ]; then
sleep 20s && echo yers && touch a

fi
if [ $SGE_TASK_ID -eq 2 ]; then
date 

fi
