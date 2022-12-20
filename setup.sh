#!/bin/bash
echo "Se crea la carpeta R_deps para instalar las librerías"
echo
mkdir -p ./R_deps
echo "Carpeta creada con éxito"
echo
echo "Llamamos a install_packages.R para cargar las librerías necesarias"
echo
Rscript --vanilla ./code/install_packages.R

echo "Librerías instaladas satisfactoriamente"

