#!/bin/bash

# This sets up the test, runs it and analyses the results.
# Any deviation from expected behaviour should be printed to stdout
# Script should return 0 on expected behaviour, 1 otherwise
# (default behaviour is to output 1 if any commands failed)


./cleanup.sh

../../sausages short.xyzclcpcs.zip test.params > stdout.output || echo "Program exited unsuccessfully"

for f in *expected
do
	diff -u $(basename $f .expected) $f
	err=$((err+$?)) # Adds return value of last command (breaks on negative returns)
done

exit $err
