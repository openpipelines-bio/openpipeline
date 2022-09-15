#!/bin/bash



echo ">>> Testing files mode"

sample="tinygex_S1_L001_I1_001"

# don't specify the mode 'files' as this is the default
./fastqc --input "cellranger_tiny_fastq/$sample.fastq.gz" --output filemode-report

echo ">> Checking whether output dir exists"
[[ ! -d filemode-report ]] && echo "Output dir could not be found!" && exit 1

echo ">> Checking if the correct files are present"
[[ ! -f filemode-report/"$sample"_fastqc.html ]] && echo "Report file missing" && exit 1
[[ ! -f filemode-report/"$sample"_fastqc.zip ]] && echo "Zip file missing" && exit 1

echo ">>> Testing dir mode"

./fastqc -m dir --input "cellranger_tiny_fastq/" --output dirmode-report

echo ">> Checking whether output dir exists"
[[ ! -d dirmode-report ]] && echo "Output dir could not be found!" && exit 1

echo ">> Checking if sufficient files are present"
# each fastq files generates one html and one zip file
nr_fastqs=`ls cellranger_tiny_fastq/*.fastq.gz | wc -l`
nr_htmlfiles=`ls dirmode-report/*.html | wc -l`
nr_zipfiles=`ls dirmode-report/*.zip | wc -l`

[[ ! $nr_fastqs == $nr_htmlfiles ]] && echo "Html files are missing" && exit 1
[[ ! $nr_fastqs == $nr_htmlfiles ]] && echo "Zip files are missing" && exit 1

# print final message
echo ">>> Test finished successfully"

# do not remove this
# as otherwise your test might exit with a different exit code
exit 0
