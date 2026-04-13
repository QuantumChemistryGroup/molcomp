#!/bin/bash

Molcomp=~/bin/Molcomp # or your path to Molcomp program
xyz=~/bin/xyz         # or your path to perl script xyz
align=~/bin/align     # or your path to perl script align
jmol=~/Downloads/jmol-14.32.10/jmol.sh # or your path to jmol.sh

input_file="$1"

while read -r file1 file2 rest; do
  name1="${file1%.xyz}"
  name2="${file2%.xyz}"
  
  echo -n "$file1 vs. $file2: "
  $Molcomp "$file1" "$file2" --print=true --connect=true --alg=std
  
  ($xyz -ellips=2 -no_hide=10 "$name1-mod.xyz"; cat "$file2") | $align -AW > "$name1-$name2.xyz"
  
  rm -f "$name1-mod.xyz"
  
  $jmol --silent --nosplash -j \
  "background white; frame ALL; set frank off; calculate HBONDS; hbonds 0.05" \
  -n -w "PNG:$name1-$name2.png" -x \
  "$name1-$name2.xyz" > /dev/null 2>&1

done < "$input_file"
