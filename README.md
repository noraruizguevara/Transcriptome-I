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
zcat SRR2922712.1.clean.fq.gz SRR2922713.1.clean.fq.gz SRR2922714.1.clean.fq.gz SRR2922715.1.clean.fq.gz SRR2922716.1.clean.fq.gz SRR2922717.1.clean.fq.gz SRR2960160.1.clean.fq.gz SRR2960161.1.clean.fq.gz SRR7003712.1.clean.fq.gz SRR7003713.1.clean.fq.gz SRR7003714.1.clean.fq.gz > all.1.clean.fq ; gzip all.1.clean.fq ; 
ls -lSh ;

zcat SRR2922712.2.clean.fq.gz SRR2922713.2.clean.fq.gz SRR2922714.2.clean.fq.gz SRR2922715.2.clean.fq.gz SRR2922716.2.clean.fq.gz SRR2922717.2.clean.fq.gz SRR2960160.2.clean.fq.gz SRR2960161.2.clean.fq.gz SRR7003712.2.clean.fq.gz SRR7003713.2.clean.fq.gz SRR7003714.2.clean.fq.gz > all.2.clean.fq ; gzip all.2.clean.fq ; 
ls -lSh ;
```
