

###### Pour créer une image singularity à partir d'une image pullée de docker image

docker pull litd/docker-cellranger:v8.0.0
# (vient de https://hub.docker.com/r/litd/docker-cellranger)

# Vérifier que l'image a été chargée
docker images

# Convertir une image docker en image singulariy
singularity build image.sif docker-daemon://litd/docker-cellranger:v8.0.0




