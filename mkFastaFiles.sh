#!/bin/bash
###############################################################################$
# Version: 0.2
# Author: Jose V. Die Ramón
# Last Change: Jose
# Date: Sun, Jul 22 09:27, 2020
#
# Take a FASTA file and make several files with N sequences per file 
#
###############################################################################$

#VARS
fastaFile=$1
nperFile=$2

printf "\n"
printf "\n"
printf "\n"

# Number of sequences in the file
nseqs=$(grep ">" $fastaFile | wc -l)
echo "The file" $fastaFile "contains" $nseqs "sequences."

# Number of required files
nfiles=$(echo "scale = 2; $nseqs / $nperFile" | bc)
echo "Number of files needed:" $nfiles

# Indexes
let counter=$nperFile+1           #1a seq del siguiente archivo
i=1                               #nombre archivos que se crean
foo=$(echo $fastaFile | cut -d '.' -f1)

start=1 #idx para 1a seq de cada archiv que se crea
firstSeq=$(head -n 1 $fastaFile | cut -d '>' -f2)
#echo -e "La 1a seq es" $firstSeq
line_firstSeq=$(grep -n $firstSeq $fastaFile | cut -d ':' -f1)
#echo -e "La linea de 1a seq es:" $line_firstSeq

#Make new directory for results
mkdir results

while [ $counter -le $nseqs ]
do
	printf "\n"
	echo -e "Making file n."$i "..."
        echo "Seq desde" $start "hasta" $counter

	# File names
        newfile=$foo$i".fasta"
        echo -e "File name:" $newfile


	#Next file
	#First seq
        target=$(grep ">" $fastaFile | head -n $counter | tail -n 1 | cut -d '>' -f2)
        echo "Next file will start with seq" $target	
        #Line for seq
        line_nextFile=$(grep -n $target $fastaFile | cut -d ':' -f1)

        #Make files
        echo -e "lineas from:" $line_firstSeq "to:" $line_nextFile
        sed -n $line_firstSeq,${line_nextFile}p $fastaFile > results/$newfile

	#Update indexes
        let start=$counter
        let counter=$counter+nperFile
        let i=$i+1
        let line_firstSeq=$line_nextFile



	#Make last file with remaining sequences
        if [ $counter -gt $nseqs ];then

          printf "\n"
                      echo "Making file n."$i "..."
                      newfile=$foo$i".fasta"
          echo -e "File name:" $newfile
                      echo -e "Contador es" $counter "pero solo hay" $nseqs
                      echo "Seq desde" $start "hasta" $nseqs
                      #echo "From line" $line_nextFile "to the end of file"

          #esta forma es la correcta pero no consigo hacerlo
                      #sed -n "$line_nextFile,${}p" $fastaFile > results/$newfile
                      jose="'640,\$p'"
                      #echo $jose
                      #sed -n ${jose} $fastaFile > results/$newfile


          #aproximacion alternativa :
                      #ultima secuencia
          lastseq=$(grep ">" $fastaFile | tail -n 1 | cut -d '>' -f2)
                      #linea de ultima secuencia
                      line_lastseq=$(grep -n $lastseq $fastaFile | cut -d ':' -f1)
          #linea final del archivo
          let file_lastLine=$line_lastseq+1

          #extract content from X to file_lastLine	
          sed -n $line_nextFile,${file_lastLine}p $fastaFile > results/$newfile


       fi

done

