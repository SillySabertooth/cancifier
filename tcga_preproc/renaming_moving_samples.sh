#it's possible to get all needed samples from TCGA https://portal.gdc.cancer.gov/
#manually by adding certain datasets in the basket
#or using manifest file via gdc tool
#we did both ways since it was not obvious in the beginning
#but after obtaining, ib both cases they should be renamed to match our data table
#on that way I used shell for some reason
#for exmple, here is a small manifest for thca datset

#mkdir ~/Tartu_Ulikooli/Intro_DS/project/data

tail -n +2 ../gdc_sample_sheet.2020-11-14.tsv | while read i; do
  dir=`echo $i | cut -f1 -d" "`
  file=`echo $i | cut -f2 -d" "`
  suf=`echo $i | cut -f2  -d" " | cut -f2,3,4 -d "."`
  name=`echo $i | cut -f10 -d" "`
#  echo $dir
#  echo $file
#  echo $suf
#  echo $name
  mv $dir/$file ~/Tartu_Ulikooli/Intro_DS/project/data/$name.$suf
done

###
#prepare files
#file=
