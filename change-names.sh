# This script is to rename the antiSMASH gbks for them to include the species and strain names, taken from the directory name.
# The argument it requires is the name of the directory with the AntiSMASH output, which must NOT contain a slash at the end.

# Usage for one AntiSMASH output directory:
#     sh change-names.sh <folder>

# Usage for multiple AntiSMASH output directory:
#     for species in <output-directory-pattern*>
#    	 do
#    		 sh change-names.sh $species
#    	 done



ls -1 "$1"/*region*gbk | while read line # enlist the gbks of all regions in the directory and start a while loop
 do
    dir=$(echo $line | cut -d'/' -f1) # save the directory name in a variable
    file=$(echo $line | cut -d'/' -f2) # save the file name in a variable
	for directory in $dir # start a for loop
    	do
   		 cd $directory # enter the directory
      	newfile=$(echo $dir-$file) # make a new variable that fuses the directory name with the file name
			 	echo "Renaming" $file " to" $newfile # print a message showing the old and new file names
			 	mv $file $newfile # rename
			 	cd .. # return to main directory before it begins again
		 done
 done
