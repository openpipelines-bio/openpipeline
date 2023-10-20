#!/bin/bash

find src -type f -exec grep --color='auto' -P -n "[\x80-\xFF]" {} +