#!/bin/bash

for directory in $(ls | grep test)
do
	echo 'Running' $directory "... "
	cd ${directory}
	./run.sh && echo "OK" || echo "Failed"
	cd ..
done

