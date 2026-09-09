# Transcriptome-I


# PASO 1: modificar headers para ser leidos por TRINITY (loop)
```r
for r1 in *_1.fastq.gz
do
prefix=$(basename $r1 _1.fastq.gz)
r2=${prefix}_2.fastq.gz
echo "-----------------------------------"
echo "CLEANING SAMPLE '${prefix}' :"
echo "-----------------------------------"
echo "working on sample '${r1}' ... "
zcat $r1 | sed -e "s/^@${prefix}./@/g" | sed -e 's/\ .*/\/1/g' > ${prefix}_f.fq ;
echo "working on sample '${r2}' ... "
zcat $r2 | sed -e "s/^@${prefix}./@/g" | sed -e 's/\ .*/\/2/g' > ${prefix}_r.fq ;
done ;

echo " " ;
echo "----------------------------------------------"
echo "NOW WE ARE COMPRESSING CLEAN FASTQ-FILES"
echo "----------------------------------------------"

gzip *.fq ;
ls -lh *.gz ;

echo "DONE!"
```

# PASO 1: modificar headers para ser leidos por TRINITY (line by line)
```r
zcat SRR7003712_1.fastq.gz | sed -e "s/^@SRR7003712./@/g" | sed -e 's/\ .*/\/1/g' > SRR7003712_f.fq ; gzip SRR7003712_f.fq ; 
zcat SRR7003712_2.fastq.gz | sed -e "s/^@SRR7003712./@/g" | sed -e 's/\ .*/\/2/g' > SRR7003712_r.fq ; gzip SRR7003712_r.fq ; 
zcat SRR7003713_1.fastq.gz | sed -e "s/^@SRR7003713./@/g" | sed -e 's/\ .*/\/1/g' > SRR7003713_f.fq ; gzip SRR7003713_f.fq ; 
zcat SRR7003713_2.fastq.gz | sed -e "s/^@SRR7003713./@/g" | sed -e 's/\ .*/\/2/g' > SRR7003713_r.fq ; gzip SRR7003713_r.fq ; 
zcat SRR7003714_1.fastq.gz | sed -e "s/^@SRR7003714./@/g" | sed -e 's/\ .*/\/1/g' > SRR7003714_f.fq ; gzip SRR7003714_f.fq ; 
zcat SRR7003714_2.fastq.gz | sed -e "s/^@SRR7003714./@/g" | sed -e 's/\ .*/\/2/g' > SRR7003714_r.fq ; gzip SRR7003714_r.fq ; 
ls -lh *.gz ;
```

# PASO 2: limpiar reads con FASTP
```r
#!/usr/bin/bash
conda activate fastp ;

RED='\e[31m'
GREEN='\e[32m'
YELLOW='\e[33m'
CYAN='\e[36m'
NC='\e[0m' # No Color

for r1 in *_f.fq.gz
do
    prefix=$(basename $r1 _f.fq.gz)
    r2=${prefix}_r.fq.gz

    echo -e "${RED}Vic: Procesando muestra ${prefix} para RNA-seq...${NC}"

    fastp -i $r1 -I $r2 \
          -o ${prefix}.1.clean.fq.gz -O ${prefix}.2.clean.fq.gz \
          --cut_mean_quality 25 \
          --cut_right \
          --qualified_quality_phred 20 \
          --unqualified_percent_limit 40 \
          --detect_adapter_for_pe \
          --trim_poly_g \
          --trim_poly_x \
          --correction \
          --n_base_limit 5 \
          -l 60 \
          --thread 8 \
          -h report.${prefix}.html \
          -j report.${prefix}.json

    echo -e "${GREEN}Vic: Muestra ${prefix} procesada exitosamente${NC}"
done ;

echo -e "${CYAN}Vic: FILTRADO FINALIZADO. Moviendo archivos...${NC}"

# Mover todo a la carpeta de reads limpios
mkdir -p 3.CLEAN.READS/ ; 
mv *.clean.fq.gz *.html *.json 3.CLEAN.READS/ 2>/dev/null || true

cd 3.CLEAN.READS/ ; 
conda deactivate ;
```

# PASO 3: Ensamblaje de transcriptomas con TRINITY
```r
#!/usr/bin/bash

conda activate trinity ; 

mkdir -p fasta/ fasta.gene_trans_map/ timing/ logs/

procesar_muestra() {
    r1=$1
    prefix=$(basename "$r1" .1.clean.fq.gz)
    r2="${prefix}.2.clean.fq.gz"
    outdir="${prefix}.trinity.out"
    
    echo "========================================="
    echo "Iniciando Trinity para: $prefix"
    echo "Fecha: $(date)"
    echo "========================================="
    
    Trinity --seqType fq \
            --left "$r1" \
            --right "$r2" \
            --max_memory 100G \
            --CPU 16 \
            --output "$outdir" \
            --no_salmon
    
    # Renombrar y mover archivos
    if [ -f "${outdir}/Trinity.fasta" ]; then
        mv "${outdir}/Trinity.fasta" "${outdir}/${prefix}.Trinity.fasta"
        mv "${outdir}/Trinity.fasta.gene_trans_map" "${outdir}/${prefix}.Trinity.fasta.gene_trans_map"
        mv "${outdir}/Trinity.timing" "${outdir}/${prefix}.Trinity.timing"
        
        mv "${outdir}/${prefix}.Trinity.fasta" fasta/
        mv "${outdir}/${prefix}.Trinity.fasta.gene_trans_map" fasta.gene_trans_map/
        mv "${outdir}/${prefix}.Trinity.timing" timing/
        
        echo "$prefix completado exitosamente"
    else
        echo "ERROR: Trinity falló para $prefix"
        exit 1
    fi
}

export -f procesar_muestra

# Procesar máximo 2 muestras en paralelo
# (cada una usa 16 núcleos para total 32 de núcleos)
ls *.1.clean.fq.gz | parallel -j 2 --bar --joblog trinity_parallel.log procesar_muestra {}

echo "DONE!"

ls -lh fasta/
```

# PASO 4: Estimacion del numero de contigs por transcriptoma ensamblado
```r

#!/usr/bin/bash

dir="4.TRINITY/fasta/"

for i in ${dir}*.fasta
do
s1=$(grep "^>" $i | wc -l)

echo "sample $i contains '$s1' contigs"
echo ""
done
```

# PASO 5: Combinar los archivos FASTQ limpios FORWARD y REVERSE para obtener el PAN-TRANSCRITOMA
```r
zcat SRR2922712.1.clean.fq.gz SRR2922713.1.clean.fq.gz SRR2922714.1.clean.fq.gz SRR2922715.1.clean.fq.gz SRR2922716.1.clean.fq.gz SRR2922717.1.clean.fq.gz SRR2960160.1.clean.fq.gz SRR2960161.1.clean.fq.gz SRR7003712.1.clean.fq.gz SRR7003713.1.clean.fq.gz SRR7003714.1.clean.fq.gz > all.1.clean.fq ;
gzip all.1.clean.fq ; 
ls -lSh ;

zcat SRR2922712.2.clean.fq.gz SRR2922713.2.clean.fq.gz SRR2922714.2.clean.fq.gz SRR2922715.2.clean.fq.gz SRR2922716.2.clean.fq.gz SRR2922717.2.clean.fq.gz SRR2960160.2.clean.fq.gz SRR2960161.2.clean.fq.gz SRR7003712.2.clean.fq.gz SRR7003713.2.clean.fq.gz SRR7003714.2.clean.fq.gz > all.2.clean.fq ;
gzip all.2.clean.fq ; 
ls -lSh ;
```

# PASO 6: Reducir la redundancia de transcritos con CD-HIT (95%)
```r
#!/usr/bin/bash

######################
## clean contigs headers ##
######################

for i in *.fasta.gz
do
p1=$(basename $i .Trinity.fasta.gz)

zcat $i | sed -e "s/^>/>${p1}_/g" | sed -e "s/\ len.*//g" | gzip > ${p1}.transcriptome.fasta.gz

p2=$(zcat $i | sed -e "s/^>/>${p1}_/g")
p3=$(zcat $i | grep "^>" | head -n 2)
p4=$(zcat $i | sed -e "s/^>/>${p1}_/g" | sed -e "s/\ len.*//g"| grep "^>" | head -n 2)
echo "INFORMATION FOR SAMPLE '${p1}'"
echo "$p3"
echo "$p4"
echo ""
done

#####################
## reduce redundance ##
#####################

conda deactivate ; 
conda activate cdhit ; 

zcat *.transcriptome.fasta.gz > all_transcriptomes_combined.fasta ;
cd-hit-est -i all_transcriptomes_combined.fasta -o metatranscriptome_reference.fasta -c 0.95 -n 10 -T 15 -M 25000 -d 0 -g 0 -r 1 ; 
ls -lSh ;
```

# PASO 7: Inferir ORFs, peptidos con TRANSDECODER 
```r
## esperar el archivo ".transdecoder.pep" ##

conda install bioconda::transdecoder
conda activate transdecoder

TransDecoder.LongOrfs -t metatranscriptome_reference.fasta ;
TransDecoder.Predict -t metatranscriptome_reference.fasta ;
```

# PASO 8: Anotar peptidos con EGGNOGMAPPER V5
 ```r
## DESCARGAR "eggnog_proteins.dmnd.gz" y "eggnog.db.gz"
wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz;
pigz -d -p 28 eggnog.db.gz

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog_proteins.dmnd.gz;
pigz -d -p 28 eggnog_proteins.dmnd.gz

wget http://eggnog5.embl.de/download/emapperdb-5.0.2/eggnog.db.gz;
pigz -t eggnog.db.gz;
pigz -d -p 28 eggnog.db.gz;
ls ;

## ANOTAR 
 emapper.py -i metatranscriptome_reference.fasta.transdecoder.pep \
           -o eggnog_annotation \
           --cpu 28 \
           -m diamond \
           --dmnd_db /mnt/d/TESIS_MACA_2026/BIOINFORMATICS/4.EGGNOG.ANNOTATION/eggnog_proteins.dmnd \
           --data_dir /mnt/d/TESIS_MACA_2026/BIOINFORMATICS/4.EGGNOG.ANNOTATION/ \
           --sensmode fast
 ```

# PASO 9: Estimar conteos por transcrito con SALMON
```r

#!/usr/bin/bash

#####################
## Estimate counts ##
#####################

salmon index -t metatranscriptome_reference.fasta -i meta_index_1 -k 31 -p 10 ; 

for r1 in *.1.clean.fq.gz
do
p1=$(basename $r1 .1.clean.fq.gz)
r2="${p1}.2.clean.fq.gz"

echo "########################################"
echo "WORKING SAMPLE '${p1}' ('$r1' and '$r2')"
echo "########################################"

salmon quant -i meta_index_1 -l A -1 "$r1" -2 "$r2" -o "quant_meta/${p1}" -p 12 --gcBias --sketch ; 
done 
echo "DONE !"

for s1 in */quant.sf
do
p1=$(echo $s1 | sed "s/\/quant.sf//g")
p2=$(grep -v "^Name" $s1 | wc -l)
echo "sample '$p1' contains '$p2' rows"
done
 ```

# PASO 10: BLAST contigs_microarray vs metatranscriptoma
 ```r
##################################
##### BLAST cinco resultados #####
##################################

# 1. Crear base de datos
makeblastdb -in metatranscriptome_reference.fasta -dbtype nucl -out meta_db ;

# 2. Ejecutar BLAST
blastn -db meta_db \
       -query contigs.array.fasta \
       -perc_identity 90 \
       -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
       -num_threads 30 \
       -max_target_seqs 5 \
       > blast_results_1.tsv

# 3. Añadir columnas de cobertura
awk -F'\t' 'BEGIN{OFS="\t"} {qcov=($4/$13)*100; scov=($4/$14)*100; print $0, qcov, scov}' blast_results_1.tsv > blast_results_with_cov_1.tsv

# 4. Añadir encabezado
sed -i '1i qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tqlen\tslen\tqcov\tscov' blast_results_with_cov_1.tsv

# 5. Ver resultados
head -n 10 blast_results_with_cov_1.tsv
 ```
