#!/bin/bash

set -e

for directory in test*
do
	echo 'Running' $directory "... "
	cd ${directory}
	./run.sh && echo "OK" || echo "Failed"
	cd ..
done

