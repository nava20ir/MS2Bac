#!/bin/bash
echo running $1
docker run -it -v "$PWD/:/fragger_projects"  fragpipebbm java -Xmx32G  -jar /fragpipe_bin/fragpipe-24.0/fragpipe-24.0/tools/MSFragger-4.4.1/MSFragger-4.4.1/MSFragger-4.4.1.jar /fragger_projects/fragger0.params  /fragger_projects/$1
