# manual_metabarcoding
Primero, hay que activar el ambiente de conda con el siguiente comando:
conda activate python2-envpython2-env

Al activarse correctamente el ambiente de conda, el nombre del ambiente aparecerá entre paréntesis antes del prompt.
Por ejemplo: (python2-env) user@nombre_computador:~$

```bash
conda activate python2-env
```
Para desactivar el ambiente de conda

```bash
conda deactivate
```

Luego, ejecuta el script de Anacapa para realizar el primer paso, que incluye la limpieza de datos y la identificación de ASVs:
Este script limpia los datos usando DADA2 y genera una tabla de ASVs.
Asegúrate de reemplazar las rutas de los archivos con las que correspondan en tu sistema.

```bash
/bin/bash /home/users/Documentos/Genom/Anacapa_db/anacapa_QC_dada2.sh -i /home/users/Documentos/Genom/dexmul_arms_al -o /home/users/Documentos/Genom/out_anacapa -d /home/users/Documentos/Genom/Anacapa_db -a nextera -t Miseq -f /home/users/Documentos/Genom/Anacapa_db/forward_primers.txt -r /home/users/Documentos/Genom/Anacapa_db/reverse_primers.txt -e /home/users/Documentos/Genom/Anacapa_db/metabarcode_loci_min_merge_length.txt -g -l
 ```

Luego, para asignar un grupo taxonómico a cada ASV identificado, debes correr el clasificador de Anacapa:
Este paso asigna la clasificación taxonómica a cada ASV basado en la base de datos proporcionada.

```bash
/bin/bash /home/users/Documentos/Genom/Anacapa_db/anacapa_classifier.sh -o  /home/users/Documentos/Genom/out_anacapa/ -d /home/users/Documentos/Genom/Anacapa_db/ -l
```
