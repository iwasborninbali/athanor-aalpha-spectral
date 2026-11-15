#!/bin/sh
set -eu
python3 interval_scan.py > interval_scan.replay.log
cmp interval_scan.log interval_scan.replay.log
python3 interval_checker.py > interval_checker.replay.log
cmp interval_checker.log interval_checker.replay.log
python3 obstruction.py > obstruction.replay.log
cmp obstruction.log obstruction.replay.log
sha256sum -c SHA256SUMS
rm interval_scan.replay.log interval_checker.replay.log obstruction.replay.log
