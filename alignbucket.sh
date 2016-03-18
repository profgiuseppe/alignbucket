#!/bin/bash

F=$1

function process  {
	echo "Buckets identified, writing fasta files..."
	python bucketize.py $F && rm buckets.list
}

if [ "$#" -lt 1 ]; then
	echo "Usage: ./alignbucket.sh <fasta file> [<delta>]"
elif [ "$#" -eq 1 ]; then
	./alignbucket -f $F && process
elif [ "$#" -eq 2 ]; then
	./alignbucket -f $F --delta $2 && process
else
	echo "Usage: ./alignbucket.sh <fasta file> [<delta>]"
fi
