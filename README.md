# ect
### Elliptic Curves Transformer


# Transformer for elliptic curves over finite fields


## Description
The program transforms points and coefficients between curves of the type Short Weierstraß, Montgomery, Edwards and Twisted Edwards.

## Requirements
python3 (with Standard Library)  
for extract_from_json.py: tkinter (e.g. \$ sudo apt install python3-tk)

## Installation and Usage
No Installation required.  
Start via: \$ python3 main.py  
Terminate at any time via: ctrl + c  
Please review https://asciinema.org/a/1OYhGARF8iyx5CDFN8IevvPBK

## Support
For support please write a comment or open an issue.

## Disclaimer
This software comes with absolutely no warranty, to the extent permitted by applicable law.  

## To Do
The method "get_zero_of_sw_function" in line 201 in curves.py utilizes a brute-force approach to find a zero of a function.  
For big numbers commonly used in modern elliptic-curve cryptography, e.g. 2 to the power of 255 - 19, this takes a very long time.  
Some research to find a faster solution didn't yield a result.  
SageMath might have a an appropriate way.  
This only affects the transformation of Short Weierstraß to Montgomery (line 214 in "trans_sw_to_m" in curves.py).  

## Project status
Development has halted.
