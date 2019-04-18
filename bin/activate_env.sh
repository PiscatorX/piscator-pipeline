#!/bin/bash -eu

basedir=`dirname ${BASH_SOURCE}`
activate_path=${basedir/bin/}python_virtualenv/bin/activate
source ${activate_path} && echo "${activate_path}"
