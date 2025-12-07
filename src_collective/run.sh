#!/bin/sh


./booksim ${1}/load_0_1 | grep -e "latency" -e "##" > ${1}/res_0_1 &
./booksim ${1}/load_0_2 | grep -e "latency" -e "##" > ${1}/res_0_2 &
./booksim ${1}/load_0_3 | grep -e "latency" -e "##" > ${1}/res_0_3 &
./booksim ${1}/load_0_4 | grep -e "latency" -e "##" > ${1}/res_0_4 &
./booksim ${1}/load_0_5 | grep -e "latency" -e "##" > ${1}/res_0_5 &
./booksim ${1}/load_0_6 | grep -e "latency" -e "##" > ${1}/res_0_6 &
./booksim ${1}/load_0_7 | grep -e "latency" -e "##" > ${1}/res_0_7 &
./booksim ${1}/load_0_8 | grep -e "latency" -e "##" > ${1}/res_0_8 &
./booksim ${1}/load_0_9 | grep -e "latency" -e "##" > ${1}/res_0_9 &
./booksim ${1}/load_1_0 | grep -e "latency" -e "##" > ${1}/res_1_0 &

