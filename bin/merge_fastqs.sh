#!/bin/bash

#INPUT-Parsing
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) input="$2"; shift ;;
        -o|--output) outputpath="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

echo "_________________"
echo "Merging start"


echo "INPUT-VALUES:"
echo "IN-Folder: $input"
echo "OUT-Folder: $outputpath"

#INPUT-Variables
in_name="$(basename -- "$input")"

echo "in_name: $in_name"

#check validity of inputs
if [ ! -d "$input" ]; then
 echo "ERROR: Path $input does not exist, please check!"
 else echo "input-path OK"
fi

if [ ! -d "$outputpath" ]; then
  echo "WARNING: $outputpath does not exist, directory will be created"
  echo "creating path $outputpath"
  mkdir -p $outputpath
  else echo "output-path OK"
fi

# list all fastq-archive-files 
gzlist=$(find "$input" -iname "*q.gz" | sort -V)
fqlist=$(find "$input" -iname "*fq" | sort -V)
fastqlist=$(find "$input" -iname "*fastq" | sort -V)

echo "All FASTQs: $gzlist"


gzcombined="merged_to_assemble.fq.gz"
# combine these fastq-files
for name in $gzlist; do
  echo "FQ.GZ: $name"
  cat "$name" >> "$gzcombined"
  gzcounter=$gzcounter+1
done
echo "merged $gzcounter FASTQ-files"
exit

findfastq=$(find . -type f -name "*_assemble.fastq.gz")


echo "#################"
echo "unzipping and NanoFilt-ering"
$(gunzip -c "$findfastq" | NanoFilt --logfile $outputpath/$in_name"_trimming.log" -q $qualityscore -l $trimlen > $outputpath/$filename"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq")
$(rm $outputpath/$filename".fastq.gz")
echo "#################"
echo "Running Flye-assembly"
# $(flye --nano-raw $outputpath/$in_name"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq" --out-dir $outputpath"/flye_assembly" --threads 4 --asm-coverage $coverage --iterations 2 --genome-size $genomesize)
$(flye --nano-raw $outputpath/$filename"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq" --out-dir $outputpath"/flye_assembly" --threads 4 --iterations 3 )
echo "#################"
echo "running medaka"
# change model after new basecalling to r941_min_sup_g507
# $(medaka_consensus -d $outputpath"/flye_assembly/assembly.fasta" -i $outputpath/$in_name"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq" -o $outputpath"/flye_medaka" -t 2 -m r941_min_hac_g507 )
echo "#################"
echo "Running minimap"

if [ -f "$outputpath""/flye_medaka/consensus.fasta" ]; then
  #$(bowtie2-build $outputpath"/flye_medaka/consensus.fasta")
  echo "Running minimap(on medaka consensus)"
  $(minimap2 -ax map-ont $outputpath"/flye_medaka/consensus.fasta" $outputpath/$filename"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq" --secondary=no | samtools view -bS | samtools sort > $outputpath/"ONT_trimmed_q_"$qualityscore"_l_"$trimlen"_to_assembly.bam"; samtools index $outputpath"/ONT_trimmed_q_"$qualityscore"_l_"$trimlen"_to_assembly.bam")
  echo "#################"
  echo "Running Prokka (on medaka consensus)"
  $(prokka $outputpath"/flye_medaka/consensus.fasta" --outdir $outputpath"/prokka_annotation" --kingdom $kingdom --cpu 4 --prefix "PROKKA" --force)
  else
    echo "Running minimap(on flye assembly)"
    $(minimap2 -ax map-ont $outputpath"/flye_assembly/assembly.fasta" $outputpath/$in_name"_trimmed_q_"$qualityscore"_l_"$trimlen".fastq" --secondary=no | samtools view -bS | samtools sort > $outputpath/"ONT_trimmed_q_"$qualityscore"_l_"$trimlen"_to_assembly.bam"; samtools index $outputpath"/ONT_trimmed_q_"$qualityscore"_l_"$trimlen"_to_assembly.bam")
    echo "#################"
    echo "Running Prokka (on flye assembly)"
    $(prokka $outputpath"/flye_assembly/assembly.fasta" --outdir $outputpath"/prokka_annotation" --kingdom $kingdom --prefix "PROKKA" --force)
fi

# pharokkainstall=$(mamba list pharokka | wc -l)
# pharokkadbs=list=$(find ./pharokka -type d -iname "pharokkadb")
# if [ "$pharokkainstall" -gt 3 ] && [ "$kingdom" == "Virus" ]; then
#   $(pharokka.py -i $outputpath"/flye_assembly/assembly.fasta" -o $outputpath"/pharokka_annotation" -d $pharokkadbs)
# fi
echo "#################"
echo "Pipeline end"
echo "_________________"
