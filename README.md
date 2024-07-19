# htslibJZ
Some simple tools using htslib

## compile for local
```sh
## Step 1: compile htslib  for other c programs (just need to do once)
cd htslib-1.20
make clean
./configure
make
cd ..

## Step 2: compile programs you need (using vcf2table as an example)
# method 1:
gcc -Wall -g -O2 -o vcf2table2 vcf2table.c -Ihtslib-1.20/htslib -Lhtslib-1.20 -lhts
# method 2:
gcc -Wall -g -O2 -o vcf2table vcf2table.c htslib-1.20/libhts.a -lbz2 -lz -lm -llzma
```

## compile for WebAssembly
```sh
# Step 1: compile htslib  for other c programs (just need to do once)
cd htslib-1.20
make clean
emconfigure ./configure CFLAGS="-s USE_ZLIB=1 -s USE_BZIP2=1" --disable-lzma
CFLAGS="-O2 -s USE_ZLIB=1 -s USE_BZIP2=1"
LDFLAGS="$EM_FLAGS -O2 -s ERROR_ON_UNDEFINED_SYMBOLS=0"
emmake make tabix CC=emcc AR=emar CFLAGS="$CFLAGS" LDFLAGS="$LDFLAGS"
cd ..

# Step 2: compile programs you need (using vcf2table as an example)
emcc $CFLAGS $LDFLAGS -o vcf2table.js vcf2table.c htslib-1.20/libhts.a -lbz2 -lz
```