#! /bin/sh

for bed in "20" "40" "60" "80" "Full"
do
for meme in "1" "178" "285" "356"
do
echo $bed.$meme
date '+TIME:%H:%M:%S'
bash -c "time python3 ./mrsoftware/mrsoftware.py -m ./tests/$meme.meme -g ./tests/GRCm38.fa ./tests/$bed.bed -o ./results/$bed$meme.txt" > ./results/$bed$meme.printed.txt
date '+TIME:%H:%M:%S'
done
done

