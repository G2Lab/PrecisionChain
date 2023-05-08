#!/bin/bash
LTIME=$(stat -c %Z /gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node/permissions.txt)
while true; do
    ATIME=$(stat -c %Z /gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node/permissions.txt)
    if [[ "$ATIME" != "$LTIME" ]]; then
        echo "RUN COMMAND"
        cat /gpfs/commons/groups/gursoy_lab/ubaymuradov/main_node/permissions.txt | bash
        LTIME=$ATIME
    fi
    sleep 2
done
