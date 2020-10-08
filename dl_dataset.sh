#!/bin/bash

# Required for large dataset that google does not analyze with an antivirus
download(){
  filename="$1"
  fileid="$2"
  curl -c ./cookie -s -L "https://drive.google.com/uc?export=download&id=${fileid}" > /dev/null
  curl -Lb ./cookie "https://drive.google.com/uc?export=download&confirm=`awk '/download/ {print $NF}' ./cookie`&id=${fileid}" -o ${filename}
  rm cookie
}

mkdir -p DATASETS
cd DATASETS;

## "small enough"
curl -L -o PAMAP2.zip 'https://drive.google.com/uc?export=download&id=0B_f2z9MW4lymdDh1aDlBSUdfbFk'
curl -L -o SoccerPos.zip 'https://drive.google.com/uc?export=download&id=0B_f2z9MW4lymaXdOMjJUa01ESUk'
curl -L -o FOG_LEG.zip 'https://drive.google.com/uc?export=download&id=0B_f2z9MW4lymblNxMmxQM0xsa2c'

## "Too large" for the "small enough" method
download MITDB.zip 0B_f2z9MW4lymME5ObHJOWkNxLXM
download REFIT.zip 0B_f2z9MW4lymZDZlZURYZktHbE0

## Unzip in parallel if possible
if ! command -v parallel &> /dev/null
then
  parallel --eta unzip -d {.} {} ::: PAMAP2.zip SoccerPos.zip FOG_LEG.zip MITDB.zip REFIT.zip
else
  unzip -d PAMAP2 PAMAP2.zip
  unzip -d SoccerPos SoccerPos.zip
  unzip -d FOG_LEG FOG_LEG.zip
  unzip -d MITDB MITDB.zip
  unzip -d REFIT REFIT.zip
fi

# The PPG dataset is missing, use alternative link
mkdir -p PPG
cd PPG
# Queries
curl -L -o Query1.txt 'https://drive.google.com/uc?export=download&id=1BOINAB7Fwe3e2JIcf-BbVP12ffjlokBa'
curl -L -o Query2.txt 'https://drive.google.com/uc?export=download&id=1Mfpyo0L0I5s788-W9qYcqBK2DDfHeBk0'
curl -L -o Query3.txt 'https://drive.google.com/uc?export=download&id=1JsQW7wyCfyZcPNGQAuDH0SkX8j1pKOnk'
curl -L -o Query4.txt 'https://drive.google.com/uc?export=download&id=1rBJZE8rKy4u9HJ1sesMY7mUIQZxA8EG4'
curl -L -o Query5.txt 'https://drive.google.com/uc?export=download&id=1_v95Oen6K6miIdMMuZtb5gOEq56ZnJCV'
# Data
download Data.txt 1Gl2xI_nERraf8gkHAuxI8uCMI5U5HbTZ

