#! /bin/bash
cp _plots/*.png .
zip -r myfile.zip ./*.png
rm *.png