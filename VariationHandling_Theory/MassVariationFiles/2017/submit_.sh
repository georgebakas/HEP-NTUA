#!/bin/bash

export X509_USER_PROXY=/afs/cern.ch/user/g/gbakas/.globus/myProxy_certificate
xrdcp $5 .
root $1 $2 $3 $4
