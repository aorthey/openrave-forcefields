#!/bin/bash
cp surfboard_com_0.txt surfboard_com_test.txt
cp surfboard_com_all_0.txt surfboard_com_all_test.txt
sed -i '$ d' surfboard_com_all_test.txt
sed -i '$ d' surfboard_com_test.txt
