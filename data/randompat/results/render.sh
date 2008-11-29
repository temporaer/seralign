#!/bin/sh

for i in *.dot ; do
	bn=`basename $i .dot`
	neato -n -Tpng:gd -o $bn.png $i
done
