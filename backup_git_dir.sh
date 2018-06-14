#!/bin/bash

rm -rf git_backup/
mkdir git_backup
cp .git/objects/ git_backup/
cp -r .git/objects/ git_backup/

