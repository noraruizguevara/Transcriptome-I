# Transcriptome-I

# PASO 0: ¿Qué motiva el desarrollo de estos comandos?
```r
El estudio de la expresión génica por RNA-SEQ requiere cuantificar las lecturas que provienen
de todos los genes transcritos durante un momento fisiológico del tejido, órgano o individuo.
El conocimiento estructural de los genes depende de contar idealmente con un genoma completo y
previamente anotado. Ante la ausencia de ello, los genes y sus iso-formas se deben inferir a 
través de sus transcritos, por ello se debe realizar un ensamblaje de-novo a partir de las
lecturas RNA-SEQ (en este caso) de diferentes muestras que corresponden a órganos de la planta.
Idealmente, y dado que los genes se expresan de manera diferencial bajo numerosas condiciones,
inferir TODOS los genes a partir de RNA-SEQ, requeriría secuenciar los ARNm de numerosas fuentes
para la misma especie (diferentes individuos, organos, tejidos, condiciones).   
```

# PASO 1: modificar headers para ser leidos por TRINITY (loop)
```r

##  Trinity tiene requerimientos para el input. En cada FASTQ file,
##  los headers descargados de NCBI necesitan ser modificados
##  para que terminen en "/1" y "/2", respectivamente.
## el siguiente loop se encarga de eso

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
# en su defecto, se pueden realizar los cambios en
# cada uno de los archivos, con las siguientes lineas

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

## Los arhivos FASTQ necesitan ser filtrados, aquellos con baja calidad seran removidos
## y se crean nuevos archivos FASTQ para luego ser ensamblados con FASTP.
## los parametros seteados en este programa tiene por objetivo :
## retener lecturas con calidades promedio de 25,
## remover los extremos 5' y 3' de baja calidad,
## retener son reads con longitudes mayores a 60 nt
## emplear 8 nucleos por operacion

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

##  Se desarrollo el siguiente codigo para ensamblar
##  dos muestras en paralelo, para ello se empleo "parallel"
## TRINITY es uno de los programas mas rigurosos para el
## ensamblaje de transcriptomas y recibe su nombre porque
## comprende 3 operaciones llamadas : "worm", "chrisalid" y "butterfly".
## el resultado mas importante es un archivo FASTA que contiene todas
## las isoformas inferidas.

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

## Es importante estimar el numero de transcritos que se han
## generado en cada una de las muestras, idealmente
## deben ser cientos de miles para Lepidium meyenii (maca)

#!/usr/bin/bash

dir="4.TRINITY/fasta/"

for i in ${dir}*.fasta
do
s1=$(grep "^>" $i | wc -l)

echo "sample $i contains '$s1' contigs"
echo ""
done
```

# PASO 5: Obtencion de un META-TRANSCRIPTOMA. Reducir la redundancia de transcritos con CD-HIT (95%)
```r

## Para obtener un meta-transcriptoma representativo y sin redundancias
## se limpian los headers de cada transcriptoma individual,
## se concatenan y finalmente se seleccionan aquellas secuencias
## (transcritos) que tienen identidades menores al 95%. El siguiente
## comando reduce la longitud de los headers de cada transcriptoma,
## genera archivos nuevos, los concatena y genera un nuevo archivo general
## que contiene secuencias con % de idetidad menores al 95%.
## de esta manera se obtiene un meta-transciptoma de referencia

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

# PASO 6: Obtencion de un PAN-TRANSCRIPTOMA. Combinar los archivos FASTQ limpios FORWARD y REVERSE.
```r

## Con el objetivo de obtener OTRO transcriptoma representativo
## de la diversidad de las 11 muestras de RNA-SEQ, empleamos otro
## procedimiento con la finalidad de obtener un pan-transcriptoma
## a partir de archivos FORWARD y REVERSE que resultan
## de la concatenacion de todos los archivos filtrados previamente.
## Las siguientes lineas de comando concatenan todos los F en un
## archivo llamado "all.1.clean.fq" que finalmente es zipeado; y 
## todos los R en un archivo llamado "all.2.clean.fq" que tambien es zipeado.
## Luego empleamos TRINITY como en el "PASO 3", pero sin necesidad de
## un "loop" porque se trata de un par de archivos para un transcriptoma unico.

zcat SRR2922712.1.clean.fq.gz SRR2922713.1.clean.fq.gz SRR2922714.1.clean.fq.gz SRR2922715.1.clean.fq.gz SRR2922716.1.clean.fq.gz SRR2922717.1.clean.fq.gz SRR2960160.1.clean.fq.gz SRR2960161.1.clean.fq.gz SRR7003712.1.clean.fq.gz SRR7003713.1.clean.fq.gz SRR7003714.1.clean.fq.gz > all.1.clean.fq ;
gzip all.1.clean.fq ; 
ls -lSh ;

zcat SRR2922712.2.clean.fq.gz SRR2922713.2.clean.fq.gz SRR2922714.2.clean.fq.gz SRR2922715.2.clean.fq.gz SRR2922716.2.clean.fq.gz SRR2922717.2.clean.fq.gz SRR2960160.2.clean.fq.gz SRR2960161.2.clean.fq.gz SRR7003712.2.clean.fq.gz SRR7003713.2.clean.fq.gz SRR7003714.2.clean.fq.gz > all.2.clean.fq ;
gzip all.2.clean.fq ; 
ls -lSh ;
```

# PASO 7: Inferir ORFs, peptidos con TRANSDECODER 
```r

## Se infieren los ORFs más largos en cada "transcriptoma
## de referencia" y se traducen para obtener secuencias peptidicas.
## En este proceso se pueden generar dos a más peptidos
## por cada transcrito. el resultado es un archivo
## de extension "*.transdecoder.pep".
## en este ejemplo se emplea solo el "meta-transcriptoma",
## pero tambien se debe emplear para el "pan-transcriptoma"

## conda install bioconda::transdecoder
## conda activate transdecoder

TransDecoder.LongOrfs -t metatranscriptome_reference.fasta ;
TransDecoder.Predict -t metatranscriptome_reference.fasta ;
```

# PASO 8: Anotar peptidos con EGGNOGMAPPER V5
 ```r

## Para anotar los peptidos se emplea EGGNOGMAPPER,
## para ello se deben descargar "eggnog_proteins.dmnd.gz" y "eggnog.db.gz"

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

## Finalmente podemos estimar las lecturas que corresponden
## a cada transcrito, es decir una tabla de conteos,
## para ello emplealos SALMON, se generan archivos llamados "quant.sf"

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

## Es importante identificar los homologos entre los contigs
## del ensayo previo (microarray) y los transcriptomas
## obtenidos a partir de los RNA-SEQ de NCBI

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

# PASOS SIGUIENTES: los siguientes pasos son el analisis de la expresion diferencial (DEGs)
 ```r
 ```
