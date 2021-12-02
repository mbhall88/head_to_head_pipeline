#!/usr/bin/env bash
if [ "$#" -ne 2 ]; then
    echo "[ERROR] Illegal number of parameters"
    echo "USAGE: UNAME=username PASSWD=password $0 <samples.xml> <submission.xml>"
    exit 2
fi

if [ -z "$PASSWD" ]; then
    echo "PASSWD is unset"
    echo "USAGE: UNAME=username PASSWD=password $0 <samples.xml> <submission.xml>"
    exit 2
fi

if [ -z "$UNAME" ]; then
    echo "UNAME is unset"
    echo "USAGE: UNAME=username PASSWD=password $0 <samples.xml> <submission.xml>"
    exit 2
fi

DST="https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"
SAMPLE_XML="$1"
SUBMISSION_XML="$2"
curl -u ${UNAME}:${PASSWD} \
    -F "SUBMISSION=@${SUBMISSION_XML}" \
    -F "SAMPLE=@${SAMPLE_XML}" \
    "$DST"
