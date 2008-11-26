#!/bin/sh

for i in ../../data/mutagenesis/common/s*.csv; do
	./hrmain --mutagenesis.in-file $i --output `basename $i .csv`.xml
done
