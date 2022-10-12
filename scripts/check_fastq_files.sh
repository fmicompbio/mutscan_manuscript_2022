## column 6 is 'pair1'
for f in `tail -n +2 $1 | cut -f 6`; do 
    if [[ ! -f "FASTQ/$f" ]]; then 
        echo "File $f doesn't exist"
        exit 1
    fi
done
## column 7 is 'pair2', but only for paired-end experiments
if [ `head -n 1 $1 | cut -f 7` == "pair2" ]; then
    echo "Paired-end data, checking pair2"
    for f in `tail -n +2 $1 | cut -f 7`; do 
        if [[ ! -f "FASTQ/$f" ]]; then 
            echo "File $f doesn't exist"
            exit 1
        fi
    done
fi
exit 0
