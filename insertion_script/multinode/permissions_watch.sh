#!/bin/bash

PERMISSIONS=/gpfs/commons/groups/gursoy_lab/aelhussein/blockchain/multichain/multinode/main_node

LTIME=$(stat -c %Z $PERMISSIONS/permissions.txt)
while true; do
    ATIME=$(stat -c %Z $PERMISSIONS/main_node/permissions.txt)
    if [[ "$ATIME" != "$LTIME" ]]; then
        echo "RUN COMMAND"
        cat $PERMISSIONS/permissions.txt | bash
        LTIME=$ATIME
    fi
    sleep 2
done
