#!/usr/bin/env bash
for val in $(ls | grep "\.[ch]p*")
do
clang-format --style=file:./.github/.clang-format ./"$val" > ./h_"$val" && mv ./h_"$val" ./"$val"
done

