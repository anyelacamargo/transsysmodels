#!/bin/sh


runseries ()
{
  lsysBasename="${1}"
  for d in 1 2 3 4 5 6 7 ; do
    ltransgl -geometry 800x600+0+0 -d $d -p 0,5,10 -r 0,90,0 -i ${lsysBasename}_d${d}.ppm ${lsysBasename}.trl 
  done
}

for wheatModel in singleshoot alwaysgreen senescence ; do
  runseries wheat_${wheatModel}
done

