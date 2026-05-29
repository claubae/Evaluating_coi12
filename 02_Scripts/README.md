# TFG_COI
En esta carpeta están tanto los scripts en R como en shell más importantes del TFG:
1) Generate_coi12.sh -> script que emplea awk para eliminar la 3era base de cada codón cuando frame = +2 (cambiar el módulo si frame es otro)
2) trimmed_fwd/rev.sh -> script cutadapt para buscar regiones ortólogas a la región Leray y eliminar los primers de dicha región 
4) clean_header.sh -> script para mantener solamente el accession y el taxID en cada header
5) extract_calanoida.sh -> script para extraer todas las secuencias pertenecientes al orden Calanoida de una base de datos.
6) Makelocaldb.sh -> script para crear una base de datos local en base a los archivos descargados desde MIDORI modificados
7) blast_sbatch.sh y noblock_blast.sh -> scripts de blastn (con y sin bloqueo taxonómico, respectivamente)
8) analyze_blast_results.R -> toma los outputs de BLASTn, los filtra y calcula la precisión taxonómica alcanzada de cada hit que haya pasado los filtros
9) pident_distribution.R -> filtra los resultados de BLASTn, calcula precisión taxonómica y transforma los valores de pident para la construcción de modeo (pident_model)
10) mafft_alignment.sh -> script para alinear las secuencias de Calanoida y poder calcular las distancias intra e interespecíficas. 