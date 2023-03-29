#!/bin/bash
julia ./calc025mm.jl&
julia ./calc2mm.jl
wait
exit 0
