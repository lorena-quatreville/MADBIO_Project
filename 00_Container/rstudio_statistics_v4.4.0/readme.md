
# rstudio

Based on rocker/tidyverse distribution (includes Rstudio): 
Added:
 - Rstudio
 - R 4.2.2
 - Misc packages for data manipulation, figures, reports
 - Seurat 4

# scDC packages found at: https://sydneybiox.github.io/scDC/

# Note : It is a good practice to put the name of the group and project into the image name while compiling
# in order to be able to associate easily the group and project name to the image

## Load a image from a gzip file:
docker load --input rstudio.tar.gz

# OR

## Build an image from a dockerfile

# IL faut être dans le répertoir où se trouve le Dockefile
cd /home/quatreville/Bureau/MADBIO_Project/02_WORKSPACE/01_MADBIO_Cameroon/01_scRNA_Seq/25_Statistical_tests/02_Container

docker build . -t madbio_rstudio_statistics

# "madbio_rstudio_statistics" is the image name
## Save
## Save the image to the project folder into the experiment 02_Container folder

docker save madbio_rstudio_statistics | gzip > /home/quatreville/Bureau/MADBIO_Project/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/25_Statistical_tests/02_Container/madbio_rstudio_statistics.tar.gz

## Run Rstudio
## Command to start a container based on the image.
## You have to change "yourPass" to a simple password (no space, no special characters).
## You can also change the mapped port number to 9090 (here 8787).

Give user details to internal script which sets user and permissions:


##### Créer un container (en le nommant différemment de l'image et avec une brèche dans le fichier "quatreville")

docker run --name rstudio_statistics -d -p 8787:8787 -v /home/quatreville/Bureau/MADBIO_Project:/home/quatreville -e PASSWORD=password -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) madbio_rstudio_statistics

# vérifier que le container a bien été créé:
docker ps -a

# vérifier qu'il tourne
docker start rstudio_statistics

# Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:8787





###### Pour recréer un container à partir d'un fichier tar contenant l'image

# Décompresser l'image:
gzip -d /home/quatreville/Bureau/MADBIO_Project/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/25_scRNASeq_A82_A86_A72_tp1_tp2/02_Container/madbio_rstudio_statistics.tar.gz


# Charger l'image
docker load < /home/quatreville/Bureau/MADBIO_Project/03_BIOINFO_PROJECT/01_MADBIO_Cameroon/01_scRNASeq/25_scRNASeq_A82_A86_A72_tp1_tp2/02_Container/madbio_rstudio_statistics.tar


# Vérifier que l'image a été chargée
docker images

##### Créer un container (en le nommant différemment de l'image et avec une brèche dans le fichier "quatreville")
docker run --name rstudio_statistics -d -p 8787:8787 -v /home/quatreville/Bureau/MADBIO_Project:/home/quatreville -e PASSWORD=password -e USER=$(whoami) -e USERID=$(id -u) -e GROUPID=$(id -g) madbio_rstudio_statistics

Then connect to the machine running docker (localhost) on mapped port (8787):
http://127.0.0.1:8787


