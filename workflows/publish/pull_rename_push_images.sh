#!/bin/bash
PROGNAME=${0##*/}
usage() {
  echo "Usage: ${PROGNAME} -t <tag> -r <target_registry> [-s <source_registry>] [-o organization] [-h] [-d]"
}

help_message() {
  cat <<- EOF
	${PROGNAME}
	Pull all docker containers from a repository, rename them and push them to a new repository.

	$(usage)

	Options:
	-t		Git tag of the verison of the images you would like to pull.
	-s		Source image repository uri.
	-r		Target container registry to push to.
	-o		Organization name.
	-h		Display this help message and exit.
	-d		Flag to enable dry-runs.
	EOF
}

SOURCE_REGISTRY="ghcr.io"
ORG="openpipelines-bio"
DRY_RUN=0
while getopts ":htr:ds" opt; do
  case ${opt} in
    h )  help_message
         exit 0;
         ;;
    t )  TAG=${OPTARG}
         ;;
    s )  SOURCE_REGISTRY=${OPTARG}
         ;;
    t )  TARGET_REGISTRY=${OPTARG}
        ;;
    o )  ORG=${OPTARG}
         ;;
    d )  DRY_RUN=1
         ;;
    \? ) echo "Invalid Option: -${OPTARG}" 1>&2
         help_message
         exit 1;
         ;;
    : )
         echo "Invalid Option: -${OPTARG} requires an argument" 1>&2
         help_message
         exit 1
         ;;
  esac
done

if [ $OPTIND -eq 1 ]; then
  echo "No options specified";
  usage
  exit 1
fi
shift $((OPTIND - 1))

# Run from root of repository
REPO_ROOT=$(git rev-parse --show-toplevel)
cd "$REPO_ROOT"

# List all folders in the "target/docker"
TARGET_FOLDER="$REPO_ROOT/target/docker"
if [[ ! -d "$TARGET_FOLDER" && -L "$TARGET_FOLDER" ]]; then
    echo "$TARGET_FOLDER does not exits, please build the docker containers first."
    exit 1;
fi

# Adopted from https://stackoverflow.com/questions/23356779/how-can-i-store-the-find-command-results-as-an-array-in-bash
# Alternative for mapfile for older bash versions.
echo "Looking for components in $TARGET_FOLDER"
components=()
tmpfile=$(mktemp)
find "$TARGET_FOLDER" -maxdepth 2 -mindepth 2 -type d -print0 > "$tmpfile"
while IFS=  read -r -d $'\0'; do
    component_name=${REPLY#"$TARGET_FOLDER/"}
    component_name_underscore=${component_name//\//_}
    components+=("$SOURCE_REGISTRY/$ORG/$component_name_underscore:$TAG")
done <"$tmpfile"
rm -f tmpfile

if [ ${#components[@]} -eq 0 ]; then
    echo "No components found in $TARGET_FOLDER, exiting!" 
    exit 1
fi

for i in "${components[@]}"
do
    printf "\t$i\n"
done

# Check if all images exist before pulling
echo "Checking if all component docker images can be found at $SOURCE_REGISTRY/$ORG."
function check_image_exists {
  docker manifest inspect "$1" > /dev/null 2> /dev/null
}
for i in "${components[@]}"
do
    check_image_exists "$i"
    exit_code=$?
    if [ $exit_code -eq 1 ]; then
        echo "Image with id $i does not exist, exiting..." && exit 1
    fi
done

# Actually pull the images
echo "Pulling images."
for i in "${components[@]}"
do
    if [ $DRY_RUN -eq 0 ]; then
        docker pull "$i" || {
            printf "Failed to pull image $i, exiting!"; 
            exit 1
        }
    else
        printf "\tDry run enabled, would try to pull $i\n"
    fi
done

# Re-tag docker containers
echo "Tagging docker containers" 
for i in "${components[@]}"
do
    if [ $DRY_RUN -eq 0 ]; then
        docker tag "$i" "${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}" || {
            echo "Failed to tag $i as ${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}, exiting...";
            exit 1
        }
    else
        printf "\tDry run enabled, would have renamed $i to ${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}\n"
    fi
done

# Pushing to new registry
echo "Pushing to registry"
for i in "${components[@]}"
do
    if [ $DRY_RUN -eq 0 ]; then
        docker image push "${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}" || {
            echo "Failed to push $i as ${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}, exiting...";
            exit 1
        }
    else
        printf "\tDry run enabled, would have pushed ${i//$SOURCE_REGISTRY/$TARGET_REGISTRY}\n"
    fi
done

echo "Finished!"