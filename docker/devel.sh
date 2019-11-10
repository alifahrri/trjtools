#!/bin/bash

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

POSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

NAME=trjtools

case $key in
    -n|--name)
    NAME="$2"
    shift # past argument
    shift # past value
    ;;
    *)    # unknown option
    POSITIONAL+=("$1") # save it in an array for later
    shift # past argument
    ;;
esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters

if [ -z "$NAME" -a "$NAME" != " " ]; then
        NAME=trjtools
fi

echo ${POSITIONAL[@]}
docker run -it -v $DIR/../:/app/trjtools ${POSITIONAL[@]} $NAME bash