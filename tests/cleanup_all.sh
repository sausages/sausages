#!/bin/bash

set -e

for directory in $(ls | grep test)
do
	echo -n 'Cleaning' $directory "... "
	cd ${directory}
	./cleanup.sh && echo "Done"
	cd ..
done

