#This script will change the word Unclassified from the ORGANISM lines by that of the respective strain code.
# Usage: sh correct_gbk.sh <gbk-file-to-edit>
file=$1   # gbk file annotated with prokka
strain=$(grep -m 1 "DEFINITION" $file |cut -d " " -f6,7) # Create a variable with the columns 6 and 7 from the DEFINITION line.

sed -i '/ORGANISM/{N;s/\n//;}' $file # Put the ORGANISM field on a single line.

sed -i "s/\s*Unclassified./ $strain/" $file # Substitute the word "Unclassified" with the value of the strain variable.
